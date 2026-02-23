# Projekt do analizy DNA

import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QFrame, QFileDialog, QMessageBox, QInputDialog
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont


class DNAAnalyzerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA Analyzer - PyQt6 GUI")
        self.setMinimumSize(1000, 700)

        # Główny widget centralny
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Główny layout pionowy
        main_layout = QVBoxLayout()
        central_widget.setLayout(main_layout)

        # =========================
        # 1. GÓRNY PANEL MENU
        # =========================
        menu_panel = QFrame()
        menu_layout = QHBoxLayout()
        menu_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        menu_panel.setLayout(menu_layout)

        menu_label = QLabel("Menu:")
        menu_label.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        menu_layout.addWidget(menu_label)

        menu_buttons = ["Plik", "Motywy", "NCBI", "Eksport", "Pomoc"]

        for name in menu_buttons:
            btn = QPushButton(name)
            btn.setFixedHeight(35)
            menu_layout.addWidget(btn)

        # =========================
        # 2. ŚRODKOWA CZĘŚĆ (LEWA + CENTRUM)
        # =========================
        middle_layout = QHBoxLayout()

        # -------- LEWY PANEL BOCZNY --------
        side_panel = QFrame()
        side_layout = QVBoxLayout()
        side_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        side_panel.setLayout(side_layout)
        side_panel.setFixedWidth(220)

        side_label = QLabel("Panel boczny")
        side_label.setFont(QFont("Arial", 11, QFont.Weight.Bold))
        side_layout.addWidget(side_label)

        side_buttons = [
            "Wczytaj plik",
            "Pobierz z NCBI",
            "Dodaj motyw",
            "Uruchom analizę",
            "Eskportuj CSV/PDF"
        ]

        self.load_button = None

        for name in side_buttons:
            btn = QPushButton(name)
            btn.setFixedHeight(40)
            side_layout.addWidget(btn)

            if name == "Wczytaj plik":
                self.load_button = btn

        if self.load_button is not None:
            self.load_button.clicked.connect(self.load_sequence)
        else:
            print("Błąd: przycisk 'Wczytaj plik' nie został znaleziony!")

        middle_layout.addWidget(side_panel)

        # -------- CENTRALNE OKNO ANALIZY --------
        center_panel = QFrame()
        center_layout = QVBoxLayout()
        center_panel.setLayout(center_layout)

        center_title = QLabel("Analizowana sekwencja DNA")
        center_title.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        center_layout.addWidget(center_title)

        # Przyciski sekcji
        section_layout = QHBoxLayout()
        section_buttons = [
            "1) Podgląd sekwencji",
            "2) Wybór motywów",
            "3) Wyniki analizy",
            "4) Wizualizacja",
            "5) Eskport"
        ]

        for name in section_buttons:
            btn = QPushButton(name)
            btn.setFixedHeight(35)
            section_layout.addWidget(btn)

        center_layout.addLayout(section_layout)

        # Okno z sekwencją DNA (placeholder)
        self.sequence_view = QTextEdit()
        self.sequence_view.setPlaceholderText("Tutaj będzie wyświetlana sekwencja DNA...")
        center_layout.addWidget(self.sequence_view)

        middle_layout.addWidget(center_panel)

        main_layout.addLayout(middle_layout)

        # =========================
        # 3. DOLNY PANEL LOGÓW
        # =========================
        log_panel = QFrame()
        log_layout = QVBoxLayout()
        log_panel.setLayout(log_layout)
        log_panel.setFixedHeight(120)

        log_label = QLabel("Logi / komunikaty")
        log_label.setFont(QFont("Arial", 11, QFont.Weight.Bold))
        log_layout.addWidget(log_label)

        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        log_layout.addWidget(self.log_output)

        main_layout.addWidget(log_panel)

    def load_sequence(self):
        options = QMessageBox(self)
        options.setWindowTitle("Wczytaj sekwencję")
        options.setText("Wybierz sposób wczytania:")

        file_button = options.addButton("Wczytaj plik FASTA", QMessageBox.ButtonRole.ActionRole)
        manual_button = options.addButton("Wpisz ręcznie", QMessageBox.ButtonRole.ActionRole)
        options.addButton("Anuluj", QMessageBox.ButtonRole.RejectRole)

        options.exec()
        clicked = options.clickedButton()

        if clicked == file_button:
            self.load_from_file()
        elif clicked == manual_button:
            self.manual_input()


    def load_from_file(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Wybierz plik FASTA",
            "",
            "FASTA Files (*.fasta *.fa *.txt)"
        )

        if not file_path:
            return

        try:
            with open(file_path, "r") as file:
                lines = file.readlines()

            sequence = ""
            for line in lines:
                if not line.startswith(">"):
                    sequence += line.strip()

            sequence = sequence.upper()

            if not self.validate_sequence(sequence):
                QMessageBox.warning(self, "Błąd", "Niepoprawne znaki DNA.")
                return

            self.sequence_view.setPlainText(sequence)
            self.log_output.append(f"Wczytano plik: {file_path}")
            self.log_output.append(f"Długość sekwencji: {len(sequence)}")

        except Exception as e:
            QMessageBox.critical(self, "Błąd", str(e))


    def manual_input(self):
        text, ok = QInputDialog.getMultiLineText(
            self,
            "Wpisz sekwencję DNA",
            "Podaj sekwencję, używając znaków A, T, C, G:"
        )

        if ok and text:
            sequence = text.replace("\n", "").replace(" ", "").upper()

            if not self.validate_sequence(sequence):
                QMessageBox.warning(self, "Błąd", "Niepoprawne litery w sekwencji DNA.")
                return

            self.sequence_view.setPlainText(sequence)
            self.log_output.append("Wprowadzono sekwencję ręcznie")
            self.log_output.append(f"Długość sekwencji: {len(sequence)}")


    def validate_sequence(self, sequence):
        return all(base in "ATCG" for base in sequence)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAAnalyzerGUI()
    window.show()
    sys.exit(app.exec())