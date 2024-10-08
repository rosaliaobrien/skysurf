import sys
from PyQt5.QtWidgets import (
    QApplication, QWidget, QGridLayout, QComboBox, QLineEdit, QLabel, QPushButton
)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import Qt
import numpy as np
import pandas as pd

class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Completeness Magnitude Calculator")
        self.layout = QGridLayout()
        self.layout.setSpacing(15)  #Increase spacing between widgets

        #Set a larger, bold font for all widgets
        font = QFont("Arial", 12)
        font.setBold(True)
        self.setFont(font)

        #Labels
        camera_label = QLabel("Camera:")
        filter_label = QLabel("Filter:")
        size_label = QLabel("Object Size (arcsec):")
        background_label = QLabel("Background (counts/sec):")
        exposure_label = QLabel("Exposure Time (sec):")

        #First ComboBox (Camera Selection)
        self.comboBox1 = QComboBox()
        self.comboBox1.addItems(["ACS", "UVIS", "IR"])
        self.comboBox1.currentIndexChanged.connect(self.update_comboBox2)

        #Second ComboBox (Filter Selection)
        self.comboBox2 = QComboBox()

        #Line edits for input
        self.object_size = QLineEdit()
        self.object_size.setPlaceholderText("Enter object size in arcseconds")

        self.background = QLineEdit()
        self.background.setPlaceholderText("Enter background value in counts/seconds")

        self.exposure_time = QLineEdit()
        self.exposure_time.setPlaceholderText("Enter exposure time in seconds")

        #Calculate button
        self.calculate_button = QPushButton("Calculate")
        self.calculate_button.clicked.connect(self.calculate_output)

        #Label for output
        self.output_label = QLabel("Output will be displayed here")
        self.output_label.setAlignment(Qt.AlignCenter)
        self.output_label.setStyleSheet("font-weight: bold; font-size: 16px; color: #FFFFFF;")

        #Add widgets to layout
        self.layout.addWidget(camera_label, 0, 0)
        self.layout.addWidget(self.comboBox1, 0, 1)

        self.layout.addWidget(filter_label, 1, 0)
        self.layout.addWidget(self.comboBox2, 1, 1)

        self.layout.addWidget(size_label, 2, 0)
        self.layout.addWidget(self.object_size, 2, 1)

        self.layout.addWidget(background_label, 3, 0)
        self.layout.addWidget(self.background, 3, 1)

        self.layout.addWidget(exposure_label, 4, 0)
        self.layout.addWidget(self.exposure_time, 4, 1)

        self.layout.addWidget(self.calculate_button, 5, 0, 1, 2)
        self.layout.addWidget(self.output_label, 6, 0, 1, 2)

        self.setLayout(self.layout)
        self.setMinimumSize(550, 400)
        self.update_comboBox2()

        #Custom Styling feel free to change
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
        """)

    def update_comboBox2(self):
        #Clear the second combobox before adding new items
        self.comboBox2.clear()

        #Update second combobox based on the first combobox camera selection
        if self.comboBox1.currentText() == "ACS":
            self.comboBox2.addItems(['F435W','F475W','F555W','F606W', 'F625W', 'F775W', 'F814W', 'F850LP'])
        elif self.comboBox1.currentText() == "UVIS":
            self.comboBox2.addItems(['F225W','F275W','F300X','F336W','F390W','F438W','F475W','F475X','F555W','F606W','F625W','F775W','F814W','F850LP'])
        elif self.comboBox1.currentText() == "IR":
            self.comboBox2.addItems(['F098M','F105W', 'F110W', 'F125W','F140W','F160W'])

    def calculate_output(self):
        #Gather input values
        try:
            camera = self.comboBox1.currentText()
            filter_number = self.comboBox2.currentText()
            wfc = '' if camera == 'ACS' else 'WFC3 '
            filter_name = f"{wfc}{camera} {filter_number}"  #Example filter name from the CSV

            background = float(self.background.text())
            exposure_time = float(self.exposure_time.text())
            size = float(self.object_size.text())

            #Load the CSV file into a DataFrame
            df = pd.read_csv(f'GUI_{camera}table.csv')
            exposure_constant = df['1" a coeff'][df.iloc[:, 0] == filter_name].values[0]
            background_constant = df['1" b coeff'][df.iloc[:, 0] == filter_name].values[0]
            avg = df['Avg'][df.iloc[:, 0] == filter_name].values[0]

            #Load elephant data for sizes
            elephant_df = pd.read_csv(f'{camera}_elephant_data.csv')
            a_elephant = elephant_df['50% a'][elephant_df.iloc[:, 0] == filter_name]
            b_elephant = elephant_df['50% b'][elephant_df.iloc[:, 0] == filter_name]
            x0_elephant = elephant_df['50% x0'][elephant_df.iloc[:, 0] == filter_name]
            c_elephant = elephant_df['50% c'][elephant_df.iloc[:, 0] == filter_name]

            if filter_number in ['F225W','F275W','F300X','F336W']: #Condition for filters with negative backgrounds. I am just ignoring their log10(background_array) component
                magnitude_constant = a_elephant + (b_elephant * (np.array(np.log10(size)) - x0_elephant - 1)) / (1 + np.exp(-c_elephant * (np.array(np.log10(size)) - x0_elephant))) - avg #np.mean(exposure_constant * np.log10(exposure_array))
            else:
                magnitude_constant = a_elephant + (b_elephant * (np.array(np.log10(size)) - x0_elephant - 1)) / (1 + np.exp(-c_elephant * (np.array(np.log10(size)) - x0_elephant))) - avg #np.mean(exposure_constant * np.log10(exposure_array) + background_constant * np.log10(background_array))

            #Calculation
            result = exposure_constant * np.log10(exposure_time) + background_constant * np.log10(background) + magnitude_constant
            result = float(result.iloc[0])

            #Display the result
            self.output_label.setText(f"Output: {result:.4f}")
        except ValueError:
            self.output_label.setText("Please enter valid numbers.")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle('Fusion')  #Set a base style

    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
