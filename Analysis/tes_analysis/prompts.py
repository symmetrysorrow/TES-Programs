import questionary


ANALYSIS_COMMANDS = [
    "Pulse Analysis",
    "Noise Analysis",
    "Temp and Optimal",
    "Scatter2D",
    "Select from Scatter",
    "ViewPulse",
    "Hist",
    "Exit",
]

KEY_CHOICES = ["Peak", "Base", "Rise", "Decay"]
MODE_CHOICES = ["Single Channel", "Two Channels"]
BIN_CHOICES = ["Auto", "Manual"]
OPTIMAL_FILTER_METHOD_CHOICES = [
    "Current (rfft/irfft + Bessel)",
    "Legacy (fft/ifft)",
]


def select_analysis_type():
    return questionary.select("Select analysis type:", choices=ANALYSIS_COMMANDS).ask()


def select_mode():
    return questionary.select("Select mode:", choices=MODE_CHOICES).ask()


def select_channel(chs):
    return questionary.select("Select Channel:", choices=chs).ask()


def select_x_key():
    return questionary.select("Select X Key:", choices=KEY_CHOICES).ask()


def select_y_key():
    return questionary.select("Select Y Key:", choices=KEY_CHOICES).ask()


def select_key():
    return questionary.select("Select Key:", choices=KEY_CHOICES).ask()


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


def input_integer(prompt):
    return questionary.text(prompt).ask()
