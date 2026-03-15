# Projekt do analizy DNA

import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QFrame, QFileDialog, QMessageBox, QInputDialog, QTabWidget, QCheckBox, QScrollArea, QStackedWidget, QTableWidget, QTableWidgetItem
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QTextCursor, QTextCharFormat, QColor
from Bio import Entrez, SeqIO
from io import StringIO
import webbrowser
import numpy as np
import pandas as pd
import random

class DNAAnalyzerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA Analyzer - PyQt6 GUI")
        self.setMinimumSize(1000, 700)
        self.motifs = []
        self.manual_sequence_counter = 1
        self.motif_colors = [
            QColor("#FFD700"),  # złoty
            QColor("#00CED1"),  # ciemny turkus
            QColor("#32CD32"),  # limonkowa zieleń
            QColor("#FF4500"),  # pomarańczowoczerwony
            QColor("#FF69B4"),  # różowy
            QColor("#1E90FF"),  # niebieski
            QColor("#8A2BE2"),  # fioletowy
            QColor("#00FA9A"),  # medium spring green
            QColor("#FF8C00"),  # dark orange
            QColor("#DA70D6")   # orchidea
]
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
            if name == "Uruchom analizę":
                btn.clicked.connect(self.run_analysis)

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

            if name == "3) Wyniki analizy":
                btn.clicked.connect(self.show_results)

        center_layout.addLayout(section_layout)

        # Okno z sekwencją DNA (placeholder)
        self.stacked_widget = QStackedWidget()
        center_layout.addWidget(self.stacked_widget)
        self.results_view = QWidget()
        self.results_layout = QVBoxLayout()
        self.results_view.setLayout(self.results_layout)

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

       # -------- Widok 3: Wyniki analizy --------
        self.analysis_view = QWidget()
        self.analysis_layout = QVBoxLayout()
        self.analysis_view.setLayout(self.analysis_layout)

        # Tabele wyników
        self.results_table = QTableWidget()
        self.segment_table = QTableWidget()

        self.analysis_layout.addWidget(QLabel("Wyniki analizy motywów:"))
        self.analysis_layout.addWidget(self.results_table)
        self.analysis_layout.addWidget(QLabel("Segmentacja sekwencji:"))
        self.analysis_layout.addWidget(self.segment_table)

        # dodanie zakładki do stacked_widget
        self.stacked_widget.addWidget(self.analysis_view)

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

            title = f"Sekwencja_{self.manual_sequence_counter}"
            self.manual_sequence_counter += 1

            self.add_sequence_tab(sequence, title=title)

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
        text_edit.setReadOnly(True)
        text_edit.setPlainText(sequence)
        layout.addWidget(text_edit)

        delete_button = QPushButton("Usuń sekwencję")
        layout.addWidget(delete_button)

        index = self.sequence_tabs.addTab(tab, title)
        self.sequence_tabs.setCurrentIndex(index)

        self.sequences.append({
            "name": title,
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

    def highlight_motifs(self, text_edit, sequence, motif_color_map):

        cursor = text_edit.textCursor()

        # reset formatowania
        cursor.select(QTextCursor.SelectionType.Document)
        cursor.setCharFormat(QTextCharFormat())

        for motif, color in motif_color_map.items():

            fmt = QTextCharFormat()
            fmt.setBackground(color)
            fmt.setForeground(QColor("black"))
            fmt.setFontWeight(QFont.Weight.Bold)

            start = 0
            motif_len = len(motif)

            while True:
                pos = sequence.find(motif, start)
                if pos == -1:
                    break

                cursor.setPosition(pos)
                cursor.movePosition(
                    QTextCursor.MoveOperation.Right,
                    QTextCursor.MoveMode.KeepAnchor,
                    motif_len
                )
                cursor.setCharFormat(fmt)

                # non-overlapping
                start = pos + motif_len

    def delete_motif_inline(self, motif_data):
        self.motifs.remove(motif_data)
        self.log_output.append(f"Usunięto motyw: {motif_data['sequence']}")
        self.show_motifs()

    def run_analysis(self):

        if not self.sequences:
            self.log_output.append("Brak sekwencji.")
            return

        active_motifs = [m["sequence"] for m in self.motifs if m["active"]]

        motif_color_map = {}
        for i, motif in enumerate(active_motifs):
            motif_color_map[motif] = self.motif_colors[i % len(self.motif_colors)]

        if not active_motifs:
            self.log_output.append("Brak aktywnych motywów.")
            return

        motif_color_map = {}
        for i, motif in enumerate(active_motifs):
            motif_color_map[motif] = self.motif_colors[i % len(self.motif_colors)]

        motif_results = []
        segment_results = []

        for seq_data in self.sequences:

            if not seq_data["checkbox"].isChecked():
                continue

            sequence = seq_data["text_edit"].toPlainText()

            self.highlight_motifs(
                seq_data["text_edit"],
                sequence,
                motif_color_map
            )
            sequence_name = seq_data["name"]

            motif_counts = self.find_motifs(sequence, active_motifs)

            for motif, count in motif_counts.items():

                motif_results.append({
                    "Sekwencja": sequence_name,
                    "Motyw": motif,
                    "Liczba": count,
                    "Długość sekwencji": len(sequence)
                })

            seg_df = self.segment_sequence(sequence)
            seg_df["Sekwencja"] = sequence_name
            seg_df = seg_df[["Sekwencja", "Początek", "Koniec", "Długość", "Zawartość GC"]]

            segment_results.append(seg_df)

        motif_df = pd.DataFrame(motif_results)

        if segment_results:
            segment_df = pd.concat(segment_results, ignore_index=True)
        else:
            segment_df = pd.DataFrame()

        self.populate_table(self.results_table, motif_df)

        if not segment_df.empty:
            self.populate_table(self.segment_table, segment_df)

        self.show_results()

        self.log_output.append("Analiza zakończona.")

    def find_motifs(self, sequence, motifs):

        results = {}

        for motif in motifs:

            count = 0
            start = 0
            motif_len = len(motif)

            while True:

                pos = sequence.find(motif, start)

                if pos == -1:
                    break

                count += 1

                start = pos + motif_len

            results[motif] = count

        return results

    def segment_sequence(self, sequence, window_size=100):

        seq_array = np.array(list(sequence))

        segments = []

        for i in range(0, len(seq_array), window_size):

            segment = seq_array[i:i+window_size]

            segment_str = "".join(segment)

            gc = (segment_str.count("G") + segment_str.count("C")) / len(segment_str)

            segments.append({
                "Początek": i,
                "Koniec": i + len(segment),
                "Długość": len(segment),
                "Zawartość GC": round(gc, 3)
            })

        df = pd.DataFrame(segments)

        return df

    def populate_table(self, table, df):

        table.clear()

        table.setRowCount(len(df))
        table.setColumnCount(len(df.columns))

        table.setHorizontalHeaderLabels(df.columns)

        for row in range(len(df)):
            for col in range(len(df.columns)):

                value = str(df.iloc[row, col])

                table.setItem(
                    row,
                    col,
                    QTableWidgetItem(value)
                )
        table.resizeColumnsToContents()

    def show_results(self):
        self.stacked_widget.setCurrentWidget(self.analysis_view)

    def display_results(self, df):

        self.results_table.clear()

        self.results_table.setRowCount(len(df))
        self.results_table.setColumnCount(len(df.columns))

        self.results_table.setHorizontalHeaderLabels(df.columns)

        for row in range(len(df)):

            for col in range(len(df.columns)):

                value = str(df.iloc[row, col])

                self.results_table.setItem(
                    row,
                    col,
                    QTableWidgetItem(value)
                )

        self.show_results()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAAnalyzerGUI()
    window.show()
    sys.exit(app.exec())