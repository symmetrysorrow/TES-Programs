#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "elmer/matc.h"

void mtc_init(FILE *, FILE *, FILE *);
char *mtc_domath(char *);

static int failures = 0;

static void expect_gate(const char *expr, int expected)
{
  int got = mtc_expression_is_temperature_only(expr);
  printf("%s\t%s\t%s\n", got == expected ? "PASS" : "FAIL",
         expected ? "accept" : "reject", expr ? expr : "<NULL>");
  if (got != expected) failures++;
}

static double eval_scalar(const char *expr)
{
  char *out = mtc_domath((char *)expr);
  char *end = NULL;
  double x;
  if (out == NULL) { failures++; return NAN; }
  x = strtod(out, &end);
  if (end == out || !isfinite(x)) {
    fprintf(stderr, "MATC expression did not return a finite scalar: %s -> %s\n", expr, out);
    failures++;
    return NAN;
  }
  return x;
}

static void expect_equivalence(const char *expr)
{
  double before = eval_scalar(expr);
  int accepted = mtc_expression_is_temperature_only(expr);
  double after = eval_scalar(expr);
  int pass = accepted && fabs(before-after) <= 1e-15 * fmax(1.0, fabs(before));
  printf("%s\teval-equivalence\t%s\t%.17g\t%.17g\n",
         pass ? "PASS" : "FAIL", expr, before, after);
  if (!pass) failures++;
}

int main(void)
{
  const char *material7 =
    "(7.854e-06)*(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))^(4.252-1)*0.4*"
    "((0.0007)-(0.0005))/(8*((1e-06)+(1.5e-05)+(1e-06))*(0.0005))";
  char before[512];
  char *last;
  double x, y;

  mtc_init(NULL, NULL, NULL);
  (void)mtc_domath("tx=0.17");
  (void)mtc_domath("guard=41");

  expect_gate("tx*tx", 1);
  expect_gate("(tx+1)^2", 1);
  expect_gate("tx-1", 1);
  expect_gate("tx-1e-12", 1);
  expect_gate("abs(tx-1)", 1);
  expect_gate("abs(tx-0.000000000001)", 1);
  expect_gate("tx/1", 1);
  expect_gate("abs(tx-1e-12)", 1);
  expect_gate("(tx+1e-12)+abs(tx-1e-12)", 1);
  expect_gate("0.5*((tx+1e-12)+abs(tx-1e-12))", 1);
  expect_gate("4.252-1", 1);
  expect_gate("(0.5*((tx+1e-12)+abs(tx-1e-12)))^(4.252-1)", 1);
  expect_gate(material7, 1);
  expect_gate("tx+ty", 0);
  expect_gate("tx+unlisted_time_state", 0);
  expect_gate("unlisted_function(tx)", 0);
  expect_gate("tx=tx+1", 0);
  expect_gate("\"text\"", 0);
  expect_gate("tx;tx+1", 0);
  expect_gate("if(tx){tx}", 0);
  expect_gate("tx+(", 0);
  expect_gate("tx 9", 0);
  expect_gate("", 0);
  expect_gate(NULL, 0);

  (void)mtc_domath("function ast_test_user(x){_ast_test_user=x*x;}");
  if (fnc_check("ast_test_user") == NULL) failures++;
  expect_gate("ast_test_user(tx)", 0);
  /* MATC resolves built-in commands before user functions and prevents a
   * user function from replacing the registered built-in. */
  (void)mtc_domath("function abs(x){_abs=x*x;}");
  if (fnc_check("abs") != NULL) {
    printf("ABS_USER_SHADOW=registered; rejecting\n");
    expect_gate("abs(tx)", 0);
  } else {
    printf("ABS_USER_SHADOW=blocked-by-MATC-builtin-namespace\n");
    expect_gate("abs(tx)", 1);
  }

  expect_equivalence("tx*tx");
  expect_equivalence("(tx+1)^2");
  expect_equivalence("abs(tx)");
  expect_equivalence(material7);

  (void)mtc_domath("guard=41");
  last = mtc_domath("guard");
  if (last == NULL) failures++;
  else { strncpy(before, last, sizeof(before)-1); before[sizeof(before)-1] = '\0'; }
  expect_gate("tx*tx", 1); /* successful parse must not replace the result buffer */
  last = mtc_domath("guard");
  if (last == NULL || strcmp(last, before) != 0) {
    fprintf(stderr, "successful parse-only check changed the prior MATC result buffer\n");
    failures++;
  }
  (void)mtc_domath("guard");
  if (fabs(eval_scalar("guard") - 41.0) > 0.0) failures++;
  last = mtc_domath("guard");
  if (last == NULL) failures++;
  else { strncpy(before, last, sizeof(before)-1); before[sizeof(before)-1] = '\0'; }
  expect_gate("tx+(", 0); /* syntax error must be fail-closed */
  last = mtc_domath("guard");
  if (last == NULL || strcmp(last, before) != 0) {
    fprintf(stderr, "parse-only check changed the prior MATC result buffer\n");
    failures++;
  }
  if (fabs(eval_scalar("guard") - 41.0) > 0.0) failures++;

  x = eval_scalar(material7);
  (void)mtc_expression_is_temperature_only(material7);
  y = eval_scalar(material7);
  if (fabs(x-y) > 1e-15 * fmax(1.0, fabs(x))) {
    fprintf(stderr, "AST check changed subsequent MATC evaluation: %.17g vs %.17g\n", x, y);
    failures++;
  }

  printf("RESULT failures=%d\n", failures);
  return failures == 0 ? 0 : 1;
}
