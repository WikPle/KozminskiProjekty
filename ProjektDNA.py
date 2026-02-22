# Projekt do analizy DNA

import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QFrame
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

        main_layout.addWidget(menu_panel)

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

        for name in side_buttons:
            btn = QPushButton(name)
            btn.setFixedHeight(40)
            side_layout.addWidget(btn)

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


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAAnalyzerGUI()
    window.show()
    sys.exit(app.exec())