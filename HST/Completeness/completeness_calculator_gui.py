import sys
from PyQt5.QtWidgets import (
    QApplication, QTabWidget, QWidget, QGridLayout, QComboBox, QLineEdit,
    QLabel, QPushButton
)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
import numpy as np
import pandas as pd

#--------------------------------------------------------------------
# Tab 1: Magnitude Calculator
#--------------------------------------------------------------------
class MagnitudeCalculatorTab(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QGridLayout()
        self.setLayout(self.layout)

        font = QFont("Arial", 12)
        font.setBold(True)
        self.setFont(font)

        camera_label = QLabel("Camera:")
        filter_label = QLabel("Filter:")
        size_label = QLabel("Object Size (arcsec):")
        background_label = QLabel("Background (counts/sec):")
        exposure_label = QLabel("Exposure Time (sec):")

        self.comboBox1 = QComboBox()
        self.comboBox1.addItems(["ACS", "UVIS", "IR"])
        self.comboBox1.currentIndexChanged.connect(self.update_comboBox2)

        self.comboBox2 = QComboBox()

        self.object_size = QLineEdit()
        self.object_size.setPlaceholderText("Enter object size in arcseconds")

        self.background = QLineEdit()
        self.background.setPlaceholderText("Enter background value in counts/sec")

        self.exposure_time = QLineEdit()
        self.exposure_time.setPlaceholderText("Enter exposure time in seconds")

        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.clicked.connect(self.calculate_output)

        self.output_label = QLabel("Output will be displayed here")
        self.output_label.setAlignment(Qt.AlignCenter)
        self.output_label.setStyleSheet("font-weight: bold; font-size: 16px; color: #FFFFFF;")

        # Place widgets in layout
        self.layout.addWidget(camera_label,     0, 0)
        self.layout.addWidget(self.comboBox1,   0, 1)
        self.layout.addWidget(filter_label,     1, 0)
        self.layout.addWidget(self.comboBox2,   1, 1)
        self.layout.addWidget(size_label,       2, 0)
        self.layout.addWidget(self.object_size, 2, 1)
        self.layout.addWidget(background_label, 3, 0)
        self.layout.addWidget(self.background,  3, 1)
        self.layout.addWidget(exposure_label,   4, 0)
        self.layout.addWidget(self.exposure_time, 4, 1)
        self.layout.addWidget(self.calculate_button, 5, 0, 1, 2)
        self.layout.addWidget(self.output_label,     6, 0, 1, 2)

        self.update_comboBox2()

    def update_comboBox2(self):
        self.comboBox2.clear()
        if self.comboBox1.currentText() == "ACS":
            self.comboBox2.addItems(['F435W','F475W','F555W','F606W',
                                     'F625W','F775W','F814W','F850LP'])
        elif self.comboBox1.currentText() == "UVIS":
            self.comboBox2.addItems(['F225W','F275W','F300X','F336W','F390W',
                                     'F438W','F475W','F475X','F555W','F606W',
                                     'F625W','F775W','F814W','F850LP'])
        elif self.comboBox1.currentText() == "IR":
            self.comboBox2.addItems(['F098M','F105W','F110W','F125W',
                                     'F140W','F160W'])

    def calculate_output(self):
        try:
            camera = self.comboBox1.currentText()
            filter_number = self.comboBox2.currentText()
            wfc = "" if camera == "ACS" else "WFC3 "
            filter_name = f"{wfc}{camera} {filter_number}"

            background = float(self.background.text())
            exposure_time = float(self.exposure_time.text())
            size = float(self.object_size.text())

            # Example CSV read
            df = pd.read_csv(f"data/{camera}table.csv")
            exposure_constant  = df['1" a coeff'][df.iloc[:, 0] == filter_name].values[0]
            background_constant= df['1" b coeff'][df.iloc[:, 0] == filter_name].values[0]
            avg                = df['Avg Exp'][df.iloc[:, 0] == filter_name].values[0]
            back                = df['Avg Back'][df.iloc[:, 0] == filter_name].values[0]

            # Elephant data (simplified)
            elephant_df = pd.read_csv(f"data/{camera}_elephant_data.csv")
            a_elephant  = elephant_df['50% a'][elephant_df.iloc[:, 0] == filter_name]
            b_elephant  = elephant_df['50% b'][elephant_df.iloc[:, 0] == filter_name]
            x0_elephant = elephant_df['50% x0'][elephant_df.iloc[:, 0] == filter_name]
            c_elephant  = elephant_df['50% c'][elephant_df.iloc[:, 0] == filter_name]

            # Example logistic formula for 50% completeness magnitude from elephant data
            magnitude_constant = (
                a_elephant 
                + (b_elephant * (np.log10(size) - x0_elephant - 1)) 
                  / (1 + np.exp(-c_elephant*(np.log10(size) - x0_elephant)))
            )

            # Combine with exposure and background terms
            result = (exposure_constant * (np.log10(exposure_time) - avg)
                      + background_constant * (np.log10(background) - back)
                      + magnitude_constant)
            result = float(result.iloc[0])
            self.output_label.setText(f"50% Completeness Mag: {result:.4f}")

        except ValueError:
            self.output_label.setText("Please enter valid numbers.")
        except FileNotFoundError:
            self.output_label.setText("CSV file not found. Check your file paths.")

#--------------------------------------------------------------------
# Tab 2: Probability (Completeness) Calculator (with SHIFT logic)
#--------------------------------------------------------------------
class ProbabilityCalculatorTab(QWidget):
    def __init__(self):
        super().__init__()
        self.layout = QGridLayout()
        self.setLayout(self.layout)

        font = QFont("Arial", 12)
        font.setBold(True)
        self.setFont(font)

        camera_label     = QLabel("Camera:")
        filter_label     = QLabel("Filter:")
        size_label       = QLabel("Object Size (arcsec):")
        background_label = QLabel("Background (counts/sec):")
        exposure_label   = QLabel("Exposure Time (sec):")
        magnitude_label  = QLabel("Object Magnitude:")

        self.comboBox1 = QComboBox()
        self.comboBox1.addItems(["ACS", "UVIS", "IR"])
        self.comboBox1.currentIndexChanged.connect(self.update_comboBox2)

        self.comboBox2 = QComboBox()

        self.object_size = QLineEdit()
        self.object_size.setPlaceholderText("Enter object size in arcseconds")

        self.background = QLineEdit()
        self.background.setPlaceholderText("Enter background in counts/sec")

        self.exposure_time = QLineEdit()
        self.exposure_time.setPlaceholderText("Enter exposure time in seconds")

        self.magnitude_input = QLineEdit()
        self.magnitude_input.setPlaceholderText("Enter object magnitude")

        self.calculate_button = QPushButton("Calculate Completeness")
        self.calculate_button.clicked.connect(self.calculate_probability)

        self.output_label = QLabel("Completeness will be displayed here")
        self.output_label.setAlignment(Qt.AlignCenter)
        self.output_label.setStyleSheet("font-weight: bold; font-size: 16px; color: #FFFFFF;")

        # Layout: Add the new Exposure Time row
        self.layout.addWidget(camera_label,     0, 0)
        self.layout.addWidget(self.comboBox1,   0, 1)
        self.layout.addWidget(filter_label,     1, 0)
        self.layout.addWidget(self.comboBox2,   1, 1)
        self.layout.addWidget(size_label,       2, 0)
        self.layout.addWidget(self.object_size, 2, 1)
        self.layout.addWidget(background_label, 3, 0)
        self.layout.addWidget(self.background,  3, 1)
        self.layout.addWidget(exposure_label,   4, 0)
        self.layout.addWidget(self.exposure_time, 4, 1)
        self.layout.addWidget(magnitude_label,  5, 0)
        self.layout.addWidget(self.magnitude_input, 5, 1)
        self.layout.addWidget(self.calculate_button, 6, 0, 1, 2)
        self.layout.addWidget(self.output_label,     7, 0, 1, 2)

        self.update_comboBox2()

    def update_comboBox2(self):
        self.comboBox2.clear()
        if self.comboBox1.currentText() == "ACS":
            self.comboBox2.addItems(['F435W','F475W','F555W','F606W',
                                     'F625W','F775W','F814W','F850LP'])
        elif self.comboBox1.currentText() == "UVIS":
            self.comboBox2.addItems(['F225W','F275W','F300X','F336W','F390W',
                                     'F438W','F475W','F475X','F555W','F606W',
                                     'F625W','F775W','F814W','F850LP'])
        elif self.comboBox1.currentText() == "IR":
            self.comboBox2.addItems(['F098M','F105W','F110W','F125W',
                                     'F140W','F160W'])

    def calculate_probability(self):
        """
        1) Compute 50% completeness from "first calculator" approach (includes exposure & background).
        2) Compute logistic-based 50%, 90%, 95% from 'elephant_data' ignoring exposure & background.
        3) Shift the 2nd calculator's values by the difference in the 50% completeness magnitudes.
        4) Fit an erfc model to the SHIFTED 50%, 90%, 95% data.
        5) Evaluate completeness at user magnitude.
        """
        try:
            from scipy.optimize import curve_fit
            from scipy.special import erfc

            # Gather user inputs
            camera = self.comboBox1.currentText()
            filter_number = self.comboBox2.currentText()
            size = float(self.object_size.text())
            background = float(self.background.text())
            exposure_time = float(self.exposure_time.text())
            mag_input = float(self.magnitude_input.text())

            #--- Read the "elephant data" CSV for the logistic size-based parameters
            elephant_df = pd.read_csv(f"data/{camera}_elephant_data.csv")
            wfc_prefix = "" if camera == "ACS" else "WFC3 "
            filter_name = f"{wfc_prefix}{camera} {filter_number}"
            row = elephant_df[elephant_df.iloc[:, 0] == filter_name]
            if row.empty:
                self.output_label.setText("No data for that Camera/Filter.")
                return

            # Elephant-based logistic parameters for 50%, 90%, 95%
            a_50, b_50, x0_50, c_50 = row["50% a"].values[0], row["50% b"].values[0], row["50% x0"].values[0], row["50% c"].values[0]
            a_90, b_90, x0_90, c_90 = row["90% a"].values[0], row["90% b"].values[0], row["90% x0"].values[0], row["90% c"].values[0]
            a_95, b_95, x0_95, c_95 = row["95% a"].values[0], row["95% b"].values[0], row["95% x0"].values[0], row["95% c"].values[0]

            def logistic_magnitude(a, b, x0, c, size_val):
                return a + (b * (np.log10(size_val) - x0 - 1)) / (
                           1 + np.exp(-c*(np.log10(size_val)-x0))
                       )

            # "elephant-based" 50%, 90%, 95% magnitudes (no exposure/background)
            mag50_elephant = logistic_magnitude(a_50, b_50, x0_50, c_50, size)
            mag90_elephant = logistic_magnitude(a_90, b_90, x0_90, c_90, size)
            mag95_elephant = logistic_magnitude(a_95, b_95, x0_95, c_95, size)

            #-----------------------------------------------------
            # 1) Compute 50% completeness from "first calculator"
            #    approach using camera table.
            #-----------------------------------------------------
            # The first calculator's formula is:
            #   mag50_firstcalc = (exposure_constant * (log10(exposure_time) - avg_exp))
            #                     + background_constant * log10(background)
            #                     + logistic_magnitude_for_50pct_size
            #   (where logistic_magnitude_for_50pct_size uses the same
            #    a_50, b_50, x0_50, c_50 from elephant data).
            df_cam = pd.read_csv(f"data/{camera}table.csv")  # same CSV used in the 1st tab
            exposure_constant  = df_cam['1" a coeff'][df_cam.iloc[:, 0] == filter_name].values[0]
            background_constant= df_cam['1" b coeff'][df_cam.iloc[:, 0] == filter_name].values[0]
            avg_exp            = df_cam['Avg Exp'][df_cam.iloc[:, 0] == filter_name].values[0]
            avg_back            = df_cam['Avg Back'][df_cam.iloc[:, 0] == filter_name].values[0]

            # logistic portion for 50% only
            mag_const_50 = mag50_elephant  # same logistic formula result above

            # Combine with exposure & background
            mag50_firstcalc = (
                exposure_constant * (np.log10(exposure_time) - avg_exp)
                + background_constant * (np.log10(background) - avg_back)
                + mag_const_50
            )

            #-------------------------------------------------------------------
            # 2) SHIFT the 50%, 90%, 95% from the elephant-based approach
            #    by Delta = (mag50_firstcalc - mag50_elephant).
            #-------------------------------------------------------------------
            delta = mag50_firstcalc - mag50_elephant
            mag50_shifted = mag50_elephant + delta
            mag90_shifted = mag90_elephant + delta
            mag95_shifted = mag95_elephant + delta

            #-------------------------------------------------------------------
            # 3) Fit an erfc model using the SHIFTED 50%, 90%, 95% data
            #-------------------------------------------------------------------
            def erfc_model(m, m0, sigma):
                # completeness(m) = 0.5 * erfc( (m - m0)/(sqrt(2)*sigma) )
                return 0.5 * erfc((m - m0)/(np.sqrt(2)*sigma))

            x_data = np.array([mag50_shifted, mag90_shifted, mag95_shifted], dtype=float)
            y_data = np.array([0.5, 0.9, 0.95], dtype=float)

            # Initial guess
            p0 = [mag50_shifted, 0.3]

            popt, pcov = curve_fit(erfc_model, x_data, y_data, p0=p0)
            m0_fit, sigma_fit = popt

            # Evaluate completeness at user magnitude
            completeness_fraction = erfc_model(mag_input, m0_fit, sigma_fit)
            completeness_percent = completeness_fraction * 100.0

            # Prepare output message
            msg = (
                f"Completeness at M={mag_input:.2f} => {completeness_percent:.2f}%\n"
                f"(Shifted) 50% Mag = {mag50_shifted:.2f}\n"
                f"(Shifted) 90% Mag = {mag90_shifted:.2f}\n"
                f"(Shifted) 95% Mag = {mag95_shifted:.2f}\n\n"
                f"First Calculator 50% Mag = {mag50_firstcalc:.2f}\n"
                f"Elephant-based 50% Mag   = {mag50_elephant:.2f}\n"
                f"Delta (shift)            = {delta:.2f}"
            )
            self.output_label.setText(msg)

        except ValueError:
            self.output_label.setText("Please enter valid numeric inputs.")
        except FileNotFoundError:
            self.output_label.setText("CSV file not found.")

#--------------------------------------------------------------------
# Main Window
#--------------------------------------------------------------------
class MainWindow(QTabWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Completeness Calculations")

        # Create the two tabs
        self.magnitude_tab = MagnitudeCalculatorTab()
        self.probability_tab = ProbabilityCalculatorTab()

        # Add them
        self.addTab(self.magnitude_tab, "Magnitude Calculation")
        self.addTab(self.probability_tab, "Probability Calculation")

        self.setMinimumSize(900, 600)  
        self.resize(900, 600)

        # Dark style
        self.setStyleSheet("""
            QWidget {
                background-color: #2E2E2E;
                color: #FFFFFF;
                font-size: 14px;
            }
            QLineEdit, QComboBox {
                background-color: #4D4D4D;
                border: 1px solid #5E5E5E;
                padding: 8px;
                font-size: 14px;
                color: #FFFFFF;
            }
            QLabel {
                font-size: 14px;
                color: #FFFFFF;
            }
            QPushButton {
                background-color: #1E90FF;
                color: #FFFFFF;
                padding: 10px;
                font-size: 16px;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #63B8FF;
            }
            QPushButton:pressed {
                background-color: #104E8B;
            }
            QComboBox QAbstractItemView {
                background-color: #4D4D4D;
                selection-background-color: #1E90FF;
                color: #FFFFFF;
            }
            QScrollBar:vertical {
                background-color: #2E2E2E;
                width: 15px;
                margin: 22px 0 22px 0;
            }
            QTabWidget::pane {
                background-color: #2E2E2E;
                border: 1px solid #444444;
            }
            QTabBar::tab {
                background-color: #4D4D4D;
                color: #FFFFFF;
                padding: 8px;
                margin-right: 2px;
            }
            QTabBar::tab:selected {
                background-color: #1E90FF;
            }
            QTabBar::tab:hover {
                background-color: #63B8FF;
            }
        """)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
