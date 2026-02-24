# Projekt do analizy DNA

import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QFrame, QFileDialog, QMessageBox, QInputDialog, QTabWidget, QCheckBox, QScrollArea, QStackedWidget
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from Bio import Entrez, SeqIO
from io import StringIO
import webbrowser


class DNAAnalyzerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA Analyzer - PyQt6 GUI")
        self.setMinimumSize(1000, 700)
        self.motifs = []
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
            "NCBI",
            "Pobierz z NCBI",
            "Dodaj motyw",
            "Uruchom analizę",
            "Eksportuj CSV/PDF",
            "Pomoc"
        ]

        self.load_button = None

        for name in side_buttons:
            btn = QPushButton(name)
            btn.setFixedHeight(40)
            side_layout.addWidget(btn)

            if name == "Wczytaj plik":
                self.load_button = btn

            if name == "Pobierz z NCBI":
                btn.clicked.connect(self.fetch_from_ncbi)

            if name == "NCBI":
                btn.clicked.connect(self.open_ncbi_browser)

            if name == "Dodaj motyw":
                btn.clicked.connect(self.add_motif)

        if self.load_button is not None:
            self.load_button.clicked.connect(self.load_sequence)
        else:
            print("Błąd: przycisk 'Wczytaj plik' nie został znaleziony!")

        middle_layout.addWidget(side_panel)

        # -------- CENTRALNE OKNO ANALIZY --------
        center_panel = QFrame()
        center_layout = QVBoxLayout()
        center_panel.setLayout(center_layout)

        center_title = QLabel("Menu analizatora DNA")
        center_title.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        center_layout.addWidget(center_title)

        # Przyciski sekcji
        section_layout = QHBoxLayout()
        section_buttons = [
            "1) Podgląd sekwencji",
            "2) Wybór motywów",
            "3) Wyniki analizy",
            "4) Wizualizacja",
        ]

        for name in section_buttons:
            btn = QPushButton(name)
            btn.setFixedHeight(35)
            section_layout.addWidget(btn)

            if name == "1) Podgląd sekwencji":
                btn.clicked.connect(self.show_sequences)

            if name == "2) Wybór motywów":
                btn.clicked.connect(self.show_motifs)

        center_layout.addLayout(section_layout)

        # Okno z sekwencją DNA (placeholder)
        self.stacked_widget = QStackedWidget()
        center_layout.addWidget(self.stacked_widget)

        # -------- Widok 1: Sekwencje --------
        self.sequence_tabs = QTabWidget()
        self.stacked_widget.addWidget(self.sequence_tabs)

        self.sequences = []

        # -------- Widok 2: Motywy --------
        self.motif_view = QWidget()
        self.motif_layout = QVBoxLayout()
        self.motif_view.setLayout(self.motif_layout)

        self.stacked_widget.addWidget(self.motif_view)

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

            self.add_sequence_tab(sequence, title="FASTA")
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

            self.add_sequence_tab(sequence, title="Manual")
            self.log_output.append("Wprowadzono sekwencję ręcznie")
            self.log_output.append(f"Długość sekwencji: {len(sequence)}")


    def validate_sequence(self, sequence):
        return all(base in "ATCG" for base in sequence)

    def open_ncbi_browser(self):
        clipboard_text = QApplication.clipboard().text().strip()

        query, ok = QInputDialog.getText(
            self,
            "Wyszukaj w NCBI Nucleotide",
            "Podaj wyszukiwany termin (np. BRCA1 human):",
            text=clipboard_text
        )

        if ok and query.strip():
            query = query.strip().replace(" ", "+")
            url = f"https://www.ncbi.nlm.nih.gov/nuccore/?term={query}"
            webbrowser.open(url)
            self.log_output.append(f"Otworzono NCBI Nucleotide dla zapytania: {query}")
        else:
            webbrowser.open("https://www.ncbi.nlm.nih.gov/nuccore/")
            self.log_output.append("Otwarto stronę główną NCBI Nucleotide.")

    def fetch_from_ncbi(self):
        clipboard = QApplication.clipboard()
        clipboard_text = clipboard.text().strip()

        accession, ok = QInputDialog.getText(
            self,
            "Pobierz z NCBI",
            "Podaj Accession ID (np. NM_001200.1):",
            text=clipboard_text
        )

        if not ok or not accession.strip():
            return

        accession = accession.strip()

        try:
            Entrez.email = "72459-CKP@kozminski.edu.pl"

            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text"
            )

            fasta_data = handle.read()
            handle.close()

            record = SeqIO.read(StringIO(fasta_data), "fasta")
            sequence = str(record.seq).upper()

            if not sequence:
                QMessageBox.warning(self, "Błąd", "Sekwencja jest pusta.")
                return

            if not self.validate_sequence(sequence):
                QMessageBox.warning(self, "Błąd", "Sekwencja zawiera niedozwolone znaki.")
                return

            self.add_sequence_tab(sequence, title=accession)

            self.log_output.append(f"Pobrano z NCBI: {accession}")
            self.log_output.append(f"Opis: {record.description}")
            self.log_output.append(f"Długość: {len(sequence)}")

        except Exception as e:
            QMessageBox.critical(self, "Błąd pobierania", str(e))

    def add_sequence_tab(self, sequence, title="Sekwencja"):
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        checkbox = QCheckBox("Aktywna sekwencja")
        checkbox.setChecked(True)
        layout.addWidget(checkbox)

        text_edit = QTextEdit()
        text_edit.setPlainText(sequence)
        layout.addWidget(text_edit)

        delete_button = QPushButton("Usuń sekwencję")
        layout.addWidget(delete_button)

        index = self.sequence_tabs.addTab(tab, title)
        self.sequence_tabs.setCurrentIndex(index)

        self.sequences.append({
            "sequence": sequence,
            "checkbox": checkbox,
            "text_edit": text_edit,
            "tab": tab
        })

        delete_button.clicked.connect(lambda: self.remove_sequence(tab))

    def show_sequences(self):
        self.stacked_widget.setCurrentWidget(self.sequence_tabs)

    def remove_sequence(self, tab):
        index = self.sequence_tabs.indexOf(tab)
        if index != -1:
            self.sequence_tabs.removeTab(index)

            self.sequences = [
                seq for seq in self.sequences
                if seq["tab"] != tab
            ]

            self.log_output.append("Usunięto sekwencję.")

    def add_motif(self):
        motif, ok = QInputDialog.getText(
            self,
            "Dodaj motyw",
            "Podaj motyw (tylko A, T, C, G):"
        )

        if not ok or not motif:
            return

        motif = motif.upper().strip()

        if not self.validate_sequence(motif):
            QMessageBox.warning(
                self,
                "Błąd",
                "Motyw może zawierać tylko litery A, T, C, G."
            )
            return

        if motif in self.motifs:
            QMessageBox.information(
                self,
                "Informacja",
                "Ten motyw już został dodany."
            )
            return

        self.motifs.append({
            "sequence": motif,
            "active": True
    })

        self.log_output.append(f"Dodano motyw: {motif}")
        self.log_output.append(f"Liczba zapisanych motywów: {len(self.motifs)}")

    def show_motifs(self):
        for i in reversed(range(self.motif_layout.count())):
            widget = self.motif_layout.itemAt(i).widget()
            if widget:
                widget.deleteLater()

        if not self.motifs:
            label = QLabel("Brak dodanych motywów.")
            self.motif_layout.addWidget(label)
        else:
            for motif_data in self.motifs:
                row = QHBoxLayout()

                checkbox = QCheckBox(motif_data["sequence"])
                checkbox.setChecked(motif_data["active"])

                checkbox.stateChanged.connect(
                    lambda state, m=motif_data: m.update({"active": state == 2})
                )

                delete_button = QPushButton("Usuń")
                delete_button.setFixedWidth(80)

                delete_button.clicked.connect(
                    lambda _, m=motif_data: self.delete_motif_inline(m)
                )

                container = QWidget()
                container.setLayout(row)

                row.addWidget(checkbox)
                row.addWidget(delete_button)

                self.motif_layout.addWidget(container)

        self.stacked_widget.setCurrentWidget(self.motif_view)

    def delete_motif_inline(self, motif_data):
        self.motifs.remove(motif_data)
        self.log_output.append(f"Usunięto motyw: {motif_data['sequence']}")
        self.show_motifs()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAAnalyzerGUI()
    window.show()
    sys.exit(app.exec())