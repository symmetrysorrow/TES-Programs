// rocALUTION PCG + AMG probe on the exported h8 CSR matrices (single GPU).
//   roc_probe A.csr Al.csr b.bin amg=rs|sa|ua|pw tilde=1 tol=3e-11 reps=3 coarse=300 acc=1
// tilde=1: AMG built on Al (preconditioning matrix), CG operator switched to A.
#include <rocalution/rocalution.hpp>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace rocalution;
typedef LocalMatrix<double> Mat;
typedef LocalVector<double> Vec;

struct Csr { int n, nnz; std::vector<int> ip, ix; std::vector<double> v; };

static Csr read_csr(const char* p)
{
    Csr m; FILE* f = fopen(p, "rb"); int h[2];
    if (!f || fread(h, sizeof(int), 2, f) != 2) { perror(p); exit(1); }
    m.n = h[0]; m.nnz = h[1]; m.ip.resize(m.n + 1); m.ix.resize(m.nnz); m.v.resize(m.nnz);
    size_t ok = fread(m.ip.data(), sizeof(int), m.n + 1, f) + fread(m.ix.data(), sizeof(int), m.nnz, f)
              + fread(m.v.data(), sizeof(double), m.nnz, f);
    (void)ok; fclose(f); return m;
}

static void load(Mat& M, const Csr& c, const char* name)
{
    M.AllocateCSR(name, c.nnz, c.n, c.n);
    M.CopyFromCSR(c.ip.data(), c.ix.data(), c.v.data());
}

static std::string kv(int argc, char** argv, const char* k, const char* d)
{
    size_t L = strlen(k);
    for (int i = 4; i < argc; i++)
        if (!strncmp(argv[i], k, L) && argv[i][L] == '=') return argv[i] + L + 1;
    return d;
}

// Preconditioner whose AMG hierarchy is built on a separate matrix (Al),
// independent of the Krylov operator (A) that CG passes in at Build().
class TildePrec : public Preconditioner<Mat, Vec, double>
{
public:
    BaseAMG<Mat, Vec, double>* amg = nullptr;
    Mat* tilde = nullptr;
    void Print(void) const override {}
    void Build(void) override { amg->SetOperator(*tilde); amg->Build(); this->build_ = true; }
    void Clear(void) override { if (amg) amg->Clear(); this->build_ = false; }
    void Solve(const Vec& rhs, Vec* x) override { amg->Solve(rhs, x); }
    void SolveZeroSol(const Vec& rhs, Vec* x) override { x->Zeros(); amg->Solve(rhs, x); }
protected:
    void PrintStart_(void) const override {}
    void PrintEnd_(void) const override {}
    void MoveToHostLocalData_(void) override {}
    void MoveToAcceleratorLocalData_(void) override {}
};

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

int main(int argc, char** argv)
{
    init_rocalution();
    info_rocalution();
    std::string amgt = kv(argc, argv, "amg", "sa");
    bool tilde = std::stoi(kv(argc, argv, "tilde", "1"));
    bool acc = std::stoi(kv(argc, argv, "acc", "1"));
    double tol = std::stod(kv(argc, argv, "tol", "3e-11"));
    int reps = std::stoi(kv(argc, argv, "reps", "3"));
    int coarse = std::stoi(kv(argc, argv, "coarse", "300"));

    Csr cA = read_csr(argv[1]), cT = read_csr(argv[2]);
    std::vector<double> hb(cA.n);
    { FILE* f = fopen(argv[3], "rb"); size_t r = fread(hb.data(), sizeof(double), cA.n, f); (void)r; fclose(f); }

    Mat A, T; Vec b, x;
    load(A, cA, "A"); load(T, cT, "Al");
    b.Allocate("b", cA.n); b.CopyFromData(hb.data());
    x.Allocate("x", cA.n);
    if (acc) { A.MoveToAccelerator(); T.MoveToAccelerator(); b.MoveToAccelerator(); x.MoveToAccelerator(); }

    for (int rep = 0; rep < reps; rep++) {
        CG<Mat, Vec, double> ls;
        Preconditioner<Mat, Vec, double>* dummy = nullptr; (void)dummy;
        BaseAMG<Mat, Vec, double>* amg;
        if (amgt == "rs") amg = new RugeStuebenAMG<Mat, Vec, double>;
        else if (amgt == "ua") amg = new UAAMG<Mat, Vec, double>;
        else if (amgt == "pw") amg = new PairwiseAMG<Mat, Vec, double>;
        else amg = new SAAMG<Mat, Vec, double>;
        amg->SetCoarsestLevel(coarse);
        amg->InitMaxIter(1);
        amg->Verbose(rep == 0 ? 2 : 0);
        x.Zeros();
        double t0 = now();
        TildePrec wrap; wrap.amg = amg; wrap.tilde = tilde ? &T : &A;
        ls.SetOperator(A);
        ls.SetPreconditioner(wrap);
        ls.Init(0.0, tol, 1e8, 2000);
        ls.Verbose(1);
        ls.Build();
        _rocalution_sync();
        double t1 = now();
        ls.Solve(b, &x);
        _rocalution_sync();
        double t2 = now();
        int its = ls.GetIterationCount();
        // true relative residual ||b - A x|| / ||b||
        Vec r; r.CloneBackend(x); r.Allocate("r", cA.n);
        A.Apply(x, &r); r.ScaleAdd(-1.0, b);
        double rel = r.Norm() / b.Norm();
        printf("amg=%s tilde=%d rep=%d build=%.3fs solve=%.3fs its=%d per_it=%.2fms true_rel=%.2e\n",
               amgt.c_str(), (int)tilde, rep, t1 - t0, t2 - t1, its, 1e3 * (t2 - t1) / (its > 0 ? its : 1), rel);
        fflush(stdout);
        ls.Clear();
        delete amg;
    }
    stop_rocalution();
    return 0;
}
