import questionary

from .analysis_utils import OPTIMAL_FILTER_METHODS


ANALYSIS_COMMANDS = [
    "Pulse Analysis",
    "Noise Analysis",
    "Temp and Optimal",
    "Compare Estimators",
    "Scatter2D",
    "Select from Scatter",
    "ViewPulse",
    "Hist",
    "Exit",
]

KEY_CHOICES = ["Peak", "Base", "Rise", "Decay"]
MODE_CHOICES = ["Single Channel", "Two Channels"]
REFINE_CHOICES = ["Finish", "Select again"]
BIN_CHOICES = ["Auto", "Manual"]
# analysis_utils が唯一の定義元。ここは表示順を決めるだけ。
OPTIMAL_FILTER_METHOD_CHOICES = list(OPTIMAL_FILTER_METHODS)


def select_analysis_type():
    return questionary.select("Select analysis type:", choices=ANALYSIS_COMMANDS).ask()


def select_mode():
    return questionary.select("Select mode:", choices=MODE_CHOICES).ask()


def select_channel(chs):
    return questionary.select("Select Channel:", choices=chs).ask()


def select_x_key(choices=None):
    return questionary.select("Select X Key:", choices=choices or KEY_CHOICES).ask()


def select_y_key(choices=None):
    return questionary.select("Select Y Key:", choices=choices or KEY_CHOICES).ask()


def select_key(choices=None):
    return questionary.select("Select Key:", choices=choices or KEY_CHOICES).ask()


def select_refine_action(count):
    return questionary.select(
        f"Selected {count} keys. Next:", choices=REFINE_CHOICES
    ).ask()


def select_two_channels(chs):
    return questionary.checkbox("Select TWO Channels:", choices=chs).ask()


def select_csv_file(files):
    return questionary.select("Select CSV File:", choices=files).ask()


def select_csv_column(columns):
    return questionary.select("Select Key:", choices=columns).ask()


def select_bin_option():
    return questionary.select("Choose bin number option:", choices=BIN_CHOICES).ask()


def select_optimal_filter_method():
    return questionary.select(
        "Select optimal filter numerical method:",
        choices=OPTIMAL_FILTER_METHOD_CHOICES,
    ).ask()


def input_eta():
    return questionary.text("eta [uA/V]:").ask()


def input_integer(prompt):
    return questionary.text(prompt).ask()
