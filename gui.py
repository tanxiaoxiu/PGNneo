import sys
import os
import threading
import time
import subprocess
from PyQt5.QtCore import QFile, QTextStream
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFileDialog,
                             QProgressBar, QTextEdit)


class UploadThread(threading.Thread):
    def __init__(self, file_path, upload_dir, progress_bar, log_text_edit):
        threading.Thread.__init__(self)
        self.file_path = file_path
        self.upload_dir = upload_dir
        self.progress_bar = progress_bar
        self.log_text_edit = log_text_edit

    def run(self):
        # Perform file upload here
        # You can replace the print statement with your own upload code
        self.log_text_edit.append(f"Uploading {self.file_path} to {self.upload_dir}...")
        for i in range(101):
            time.sleep(0.05)  # Simulate file upload
            self.progress_bar.setValue(i)
        self.log_text_edit.append(f"{self.file_path} has been uploaded to {self.upload_dir}.")


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setGeometry(300, 100, 1200, 800)
        # Initialize the UI components
        title_style = '''
            QLabel {
                font-size: 30px;
                font-weight: bold;
                color: #333;
                padding: 10px;
                qproperty-alignment: AlignCenter;
            }
        '''
        self.title_labels=[]
        self.file_labels = []
        self.file_buttons = []
        self.title=QLabel("Welcome to PGNneo")
        self.title.setStyleSheet(title_style)
        self.upload_dir_label = QLabel("Upload Directory:")
        self.upload_dir_button = QPushButton("Choose Directory")
        self.upload_dir_button.clicked.connect(self.choose_directory)
        self.upload_button = QPushButton("Start Analysis")
        self.upload_button.clicked.connect(self.upload_files)
        #self.progress_bars = []
        self.log_text_edit = QTextEdit()

        # Set the layout
        main_layout = QVBoxLayout()
        file_layouts = []
        main_layout.addWidget(self.title)
        main_layout.addWidget(self.upload_dir_label)
        main_layout.addWidget(self.upload_dir_button)
        j=['Control_R1','Control_R2','Case_R1','Case_R2','MassSpectrumData']
        for i in range(5):
            file_layouts.append(QHBoxLayout())
            self.file_labels.append(QLabel(f"{j[i]}:"))
            self.file_buttons.append(QPushButton("Choose File"))
            self.file_buttons[-1].clicked.connect(self.choose_file)
            #self.progress_bars.append(QProgressBar())
            #self.progress_bars[-1].setValue(0)
            file_layouts[-1].addWidget(self.file_labels[-1])
            file_layouts[-1].addWidget(self.file_buttons[-1])
            #file_layouts[-1].addWidget(self.progress_bars[-1])
            main_layout.addLayout(file_layouts[-1])

        main_layout.addWidget(self.upload_button)
        main_layout.addWidget(QLabel("Log:"))
        main_layout.addWidget(self.log_text_edit)
        self.setLayout(main_layout)

        # Set the upload directory to the current working directory
        self.upload_dir = os.getcwd()

    def choose_file(self):
        # Open a file dialog to select a file
        file_path, _ = QFileDialog.getOpenFileName(self, "Choose File")
        sender = self.sender()
        index = self.file_buttons.index(sender)
        self.file_labels[index].setText(file_path)
        self.log_text_edit.append(f"Your select file is {file_path}")

    def choose_directory(self):
        # Open a file dialog to select a directory
        self.upload_dir = QFileDialog.getExistingDirectory(self, "Choose Directory")
        self.upload_dir_label.setText(f"PGNneo Directory: {self.upload_dir}")
        os.chdir(self.upload_dir)  ##修改了当前路径

    def upload_files(self):
        # Upload the selected files
        self.log_text_edit.append(f"Current Wokrdir: {os.getcwd()}")
        for i in range(5):
            file_path = self.file_labels[i].text()
            #if file_path != "":
                #print(file_path)
                #print(os.getcwd())
                #print('model %d success'% i)
                #self.log_text_edit.append(f"{file_path} has been uploaded.")
                #self.log_text_edit.append(f"{file_path} has been uploaded.")
                #self.log_text_edit.append(f"Model {i} Success.")
                #thread = UploadThread(file_path, self.upload_dir, self.progress_bars[i], self.log_text_edit)
                #thread.start()
        self.log_text_edit.append(f"{file_path}")
        self.log_text_edit.append(f"Strat analysis!")
        self.log_text_edit.append(f">>>>>Model1")
        control1=self.file_labels[0].text()
        control2=self.file_labels[1].text()
        case1=self.file_labels[2].text()
        case2=self.file_labels[3].text()
        try:
            try:
                subprocess.run(['mkdir -p test'])
                subprocess.run(['cp',control1,control2,case1,case2,'./test/'])
            except Exception as e:
                self.log_text_edit.append(f"Out of Memory")
            self.log_text_edit.append(f"python model1_rnaseq_mutation_hla.py {os.path.basename(control1)} {os.path.basename(control2)} {os.path.basename(case1)} {os.path.basename(case2)}")
            try:
                result =subprocess.run(["python model1_rnaseq_mutation_hla.py", os.path.basename(control1),os.path.basename(control2),os.path.basename(case1),os.path.basename(case2)], capture_output=True, text=True).stdout
                self.log_text_edit.append(f"{result}")
            except Exception as e:
                self.log_text_edit.append(f"Fail,Please check your ENV!")
            #subprocess.run(["ls"], capture_output=True, text=True).stdout
            self.log_text_edit.append(f">>>>>Model2")
            self.log_text_edit.append(f"python model2_mutated_peptides.py {os.path.basename(control1)} {os.path.basename(case1)}")
            try:
                result =subprocess.run(["python model2_mutated_peptides.py", os.path.basename(control1),os.path.basename(case1)], capture_output=True, text=True).stdout
                self.log_text_edit.append(f"{result}")
            except Exception as e:
                self.log_text_edit.append(f"Fail,Please check your ENV!")
            self.log_text_edit.append(f">>>>>Model3")
            self.log_text_edit.append(f"python model3_MS_filtration.py")
            try:
                result =subprocess.run(["python model2_mutated_peptides.py"], capture_output=True, text=True).stdout
                self.log_text_edit.append(f"{result}")
            except Exception as e:
                self.log_text_edit.append(f"Fail,Please check your ENV!")
            self.log_text_edit.append(f">>>>>Model4")
            self.log_text_edit.append(f"python model4_neoantigen_prediction_filtration.py")
            try:
                result = subprocess.run(["python model4_neoantigen_prediction_filtration.py"], capture_output=True, text=True).stdout
                self.log_text_edit.append(f"{result}")
            except Exception as e:
                self.log_text_edit.append(f"Fail,Please check your ENV!")
            self.log_text_edit.append(f"Your results were generated in current wokrdir!")
        except Exception as e:
            self.log_text_edit.append(f"Fail,Please check your ENV!")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    file = QFile("style.qss")
    file.open(QFile.ReadOnly | QFile.Text)
    stream = QTextStream(file)
    app.setStyleSheet(stream.readAll())
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())

