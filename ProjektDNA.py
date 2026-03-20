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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages

class DNAAnalyzerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DNA Analyzer - PyQt6 GUI")
        self.setMinimumSize(1000, 700)
        self.motifs = []
        self.wykresy = Wykresy(self)
        self.report_exporter = RaportExporter(self)
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

            if name == "Eksportuj CSV/PDF":
                btn.clicked.connect(self.show_export_options)

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

            if name == "4) Wizualizacja":
                btn.clicked.connect(self.show_visualization)

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

        # zakładki wyników
        self.analysis_tabs = QTabWidget()

        # --- TAB 1 MOTIFS ---
        self.motif_results_tab = QWidget()
        motif_layout = QVBoxLayout()
        self.motif_results_tab.setLayout(motif_layout)

        self.results_table = QTableWidget()
        motif_layout.addWidget(self.results_table)
        self.analysis_tabs.addTab(self.motif_results_tab, "Motywy")

        # --- TAB 2 SUMMARY ---
        self.summary_tab = QWidget()
        self.summary_layout = QVBoxLayout()
        self.summary_tab.setLayout(self.summary_layout)

        self.summary_table = QTableWidget()
        self.summary_layout.addWidget(self.summary_table)

        self.analysis_tabs.addTab(self.summary_tab, "Podsumowanie motywów")

        # --- TAB 3 SEGMENTS ---
        self.segment_results_tab = QWidget()
        segment_layout = QVBoxLayout()
        self.segment_results_tab.setLayout(segment_layout)

        self.segment_table = QTableWidget()
        segment_layout.addWidget(self.segment_table)
        self.analysis_tabs.addTab(self.segment_results_tab, "Segmentacja GC")


        self.analysis_layout.addWidget(self.analysis_tabs)

        self.summary_label = QLabel("")
        self.summary_label.setWordWrap(True)
        self.summary_label.setStyleSheet("padding: 10px;")
        self.analysis_layout.addWidget(self.summary_label)

        self.stacked_widget.addWidget(self.analysis_view)

        # -------- Widok 4: Wizualizacja --------
        self.visualization_view = QWidget()
        self.visual_layout = QVBoxLayout()
        self.visualization_view.setLayout(self.visual_layout)

        # zakładki wizualizacji
        self.visual_tabs = QTabWidget()

        # --- TAB 1 BAR CHART ---
        self.bar_tab = QWidget()
        bar_layout = QVBoxLayout()
        self.bar_tab.setLayout(bar_layout)

        self.bar_canvas = FigureCanvas(plt.figure())
        bar_layout.addWidget(self.bar_canvas)

        # --- TAB 2 HEATMAP ---
        self.heatmap_tab = QWidget()
        heatmap_layout = QVBoxLayout()
        self.heatmap_tab.setLayout(heatmap_layout)

        self.heatmap_canvas = FigureCanvas(plt.figure())
        heatmap_layout.addWidget(self.heatmap_canvas)

        # --- TAB 3 MOTIF POSITIONS ---
        self.position_tab = QWidget()
        position_layout = QVBoxLayout()
        self.position_tab.setLayout(position_layout)

        self.position_canvas = FigureCanvas(plt.figure())
        position_layout.addWidget(self.position_canvas)

        # dodanie zakładek
        self.visual_tabs.addTab(self.bar_tab, "Wykres motywów")
        self.visual_tabs.addTab(self.heatmap_tab, "Heatmapa")
        self.visual_tabs.addTab(self.position_tab, "Pozycje motywów")

        self.visual_layout.addWidget(self.visual_tabs)

        self.stacked_widget.addWidget(self.visualization_view)

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

            sequence = seq_data["sequence"]

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

        self.motif_df = pd.DataFrame(motif_results)
        motif_df = self.motif_df

        pivot = self.motif_df.pivot_table(
            index="Motyw",
            columns="Sekwencja",
            values="Liczba",
            fill_value=0
        )

        self.pivot_df = pivot.reset_index()

        if segment_results:
            segment_df = pd.concat(segment_results, ignore_index=True)
        else:
            segment_df = pd.DataFrame()

        self.segment_df = segment_df

        self.populate_table(self.results_table, self.pivot_df)

        summary_df = self.analyze_motifs_summary(self.pivot_df)
        self.populate_table(self.summary_table, summary_df)

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

    def analyze_motifs_summary(self, pivot_df):

        df = pivot_df.set_index("Motyw")

        results = []

        # 🔹 wspólne
        common = df[(df > 0).all(axis=1)].index.tolist()

        results.append({
            "Typ": "Wspólne",
            "Sekwencja": "Wszystkie",
            "Liczba": len(common),
            "Motywy": ", ".join(common) if common else "Brak"
        })

        # 🔹 unikalne
        for col in df.columns:
            unique = df[
                (df[col] > 0) &
                (df.drop(columns=[col]) == 0).all(axis=1)
            ].index.tolist()

            results.append({
                "Typ": "Unikalne",
                "Sekwencja": col,
                "Liczba": len(unique),
                "Motywy": ", ".join(unique) if unique else "Brak"
            })

        return pd.DataFrame(results)

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

    def show_visualization(self):
        if not hasattr(self, "motif_df"):
            self.log_output.append("Najpierw uruchom analizę.")
            return

        self.wykresy.plot_bar_chart()
        self.wykresy.plot_heatmap_visualization()
        self.wykresy.plot_motif_positions_visualization()

        self.stacked_widget.setCurrentWidget(self.visualization_view)

    def show_export_options(self):
        options = QMessageBox(self)
        options.setWindowTitle("Eksport danych")
        options.setText("Wybierz opcję eksportu:")

        csv_button = options.addButton("Eksportuj CSV", QMessageBox.ButtonRole.ActionRole)
        pdf_button = options.addButton("Eksportuj PDF", QMessageBox.ButtonRole.ActionRole)
        options.addButton("Anuluj", QMessageBox.ButtonRole.RejectRole)

        options.exec()
        clicked = options.clickedButton()

        if clicked == csv_button:
            self.report_exporter.export_csv()
        elif clicked == pdf_button:
            self.report_exporter.export_pdf()

class Wykresy:
    def __init__(self, analyzer_gui):
        self.gui = analyzer_gui  # referencja do głównej klasy GUI

    def plot_bar_chart(self):
        fig = self.gui.bar_canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)

        pivot = self.gui.motif_df.pivot_table(
            index="Sekwencja",
            columns="Motyw",
            values="Liczba",
            fill_value=0
        )

        pivot.plot(kind="bar", ax=ax)
        ax.set_title("Liczba motywów w sekwencjach")
        ax.set_ylabel("Liczba wystąpień")
        ax.set_xlabel("Sekwencja")

        fig.tight_layout()
        self.gui.bar_canvas.draw()

    def plot_heatmap_visualization(self):
        active_sequences = [seq for seq in self.gui.sequences if seq["checkbox"].isChecked()]
        motifs = [m["sequence"] for m in self.gui.motifs if m["active"]]

        if not active_sequences or not motifs:
            return

        fig = self.gui.heatmap_canvas.figure
        fig.clear()
        n = len(active_sequences)
        fig.set_size_inches(10, max(4, 3 * n))

        for idx, seq_data in enumerate(active_sequences):
            sequence = seq_data["sequence"]
            name = seq_data["name"]
            window = 100
            matrix = []

            for motif in motifs:
                row = []
                for i in range(0, len(sequence), window):
                    segment = sequence[i:i+window]
                    row.append(segment.count(motif))
                matrix.append(row)
            matrix = np.array(matrix)

            ax = fig.add_subplot(n, 1, idx+1)
            im = ax.imshow(matrix, cmap="viridis", aspect="auto")

            ax.set_title(f"Heatmapa motywów – {name}", fontsize=10)
            ax.set_xlabel("Segment")
            ax.set_ylabel("Motyw")
            ax.set_yticks(range(len(motifs)))
            ax.set_yticklabels(motifs, fontsize=9)
            ax.set_xticks(range(matrix.shape[1]))
            ax.set_xticklabels(range(1, matrix.shape[1]+1), fontsize=9)

        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        cbar_ax.set_title("Liczba", fontsize=9)
        fig.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.1, hspace=0.5)
        self.gui.heatmap_canvas.draw()

    def plot_motif_positions_visualization(self):
        active_sequences = [seq for seq in self.gui.sequences if seq["checkbox"].isChecked()]
        motifs = [m["sequence"] for m in self.gui.motifs if m["active"]]

        if not active_sequences or not motifs:
            return

        fig = self.gui.position_canvas.figure
        fig.clear()
        n = len(active_sequences)
        fig.set_size_inches(10, max(4, 3 * n))

        for idx, seq_data in enumerate(active_sequences):
            sequence = seq_data["sequence"]
            name = seq_data["name"]
            ax = fig.add_subplot(n, 1, idx+1)

            y = 0
            for motif in motifs:
                positions = []
                start = 0
                while True:
                    pos = sequence.find(motif, start)
                    if pos == -1:
                        break
                    positions.append(pos)
                    start = pos + 1

                ax.scatter(positions, [y]*len(positions), label=motif)
                y += 1

            ax.set_title(f"Pozycje motywów – {name}")
            ax.set_xlabel("Pozycja w sekwencji")
            ax.set_yticks(range(len(motifs)))
            ax.set_yticklabels(motifs)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

        fig.subplots_adjust(left=0.07, right=0.8, top=0.95, bottom=0.05, hspace=0.5)
        self.gui.position_canvas.draw()

class RaportExporter:
    def __init__(self, gui):
        self.gui = gui

    def export_csv(self):
        if not hasattr(self.gui, "motif_df") or self.gui.motif_df.empty:
            QMessageBox.warning(self.gui, "Błąd", "Brak danych do eksportu.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self.gui,
            "Zapisz dane do CSV",
            "",
            "CSV Files (*.csv)"
        )

        if not file_path:
            return

        try:
            # 🔹 Zapis podsumowania statystyk
            summary_df = self.summarize_statistics()
            summary_df.to_csv(file_path, index=False)

            # 🔹 Zapis pivot_df (liczba motywów na sekwencję)
            pivot_file = file_path.replace(".csv", "_motywy.csv")
            self.gui.pivot_df.to_csv(pivot_file, index=False)

            # 🔹 Zapis segment_df (segmentacja i GC)
            if hasattr(self.gui, "segment_df") and not self.gui.segment_df.empty:
                segment_file = file_path.replace(".csv", "_segmenty.csv")
                self.gui.segment_df.to_csv(segment_file, index=False)

            QMessageBox.information(
                self.gui,
                "Eksport CSV",
                f"Dane zapisano do:\n{file_path}\n{pivot_file}\n{segment_file if hasattr(self.gui, 'segment_df') else ''}"
            )
            self.gui.log_output.append(f"Wyeksportowano dane CSV: {file_path}, pivot: {pivot_file}")

        except Exception as e:
            QMessageBox.critical(self.gui, "Błąd", f"Nie udało się zapisać CSV:\n{str(e)}")

    def export_pdf(self):
        if not hasattr(self.gui, "motif_df") or self.gui.motif_df.empty:
            QMessageBox.warning(self.gui, "Błąd", "Brak danych do raportu.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self.gui,
            "Zapisz raport do PDF",
            "",
            "PDF Files (*.pdf)"
        )

        if not file_path:
            return

        try:
            with PdfPages(file_path) as pdf:

                fig, ax = plt.subplots(figsize=(8, 6))
                ax.axis('off')
                summary_df = self.summarize_statistics()

                if summary_df.empty:
                    raise ValueError("Brak danych do podsumowania.")

                rows, cols = summary_df.shape
                fig.set_size_inches(min(20, max(8, cols * 1.2)),
                                    min(20, max(6, rows * 0.6)))

                table = ax.table(
                    cellText=summary_df.values,
                    colLabels=summary_df.columns,
                    cellLoc='center',
                    colLoc='center',
                    loc='center'
                )

                fontsize = max(6, min(12, 20 / max(rows, cols)))
                table.auto_set_font_size(False)
                table.set_fontsize(fontsize)
                table.scale(min(1.5, 20 / cols), min(1.5, 20 / rows))

                plt.tight_layout()
                ax.set_title("Podsumowanie statystyk")
                pdf.savefig(fig)
                plt.close(fig)

                fig, ax = plt.subplots(figsize=(10, 6))
                ax.axis('off')

                rows, cols = self.gui.pivot_df.shape
                fig.set_size_inches(min(20, max(8, cols * 1.2)),
                                    min(20, max(6, rows * 0.4)))

                table = ax.table(
                    cellText=self.gui.pivot_df.values,
                    colLabels=self.gui.pivot_df.columns,
                    cellLoc='center',
                    colLoc='center',
                    loc='center'
                )

                fontsize = max(6, min(12, 20 / max(rows, cols)))
                table.auto_set_font_size(False)
                table.set_fontsize(fontsize)
                table.scale(min(1.5, 20 / cols), min(1.5, 20 / rows))

                plt.tight_layout()
                ax.set_title("Liczba motywów w sekwencjach")
                pdf.savefig(fig)
                plt.close(fig)

                if hasattr(self.gui, "segment_df") and not self.gui.segment_df.empty:
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.axis('off')

                    rows, cols = self.gui.segment_df.shape
                    fig.set_size_inches(min(20, max(8, cols * 1.2)),
                                        min(20, max(6, rows * 0.4)))

                    table = ax.table(
                        cellText=self.gui.segment_df.values,
                        colLabels=self.gui.segment_df.columns,
                        cellLoc='center',
                        colLoc='center',
                        loc='center'
                    )

                    fontsize = max(6, min(12, 20 / max(rows, cols)))
                    table.auto_set_font_size(False)
                    table.set_fontsize(fontsize)
                    table.scale(min(1.5, 20 / cols), min(1.5, 20 / rows))

                    plt.tight_layout()
                    ax.set_title("Segmentacja sekwencji i zawartość GC", pad=20)
                    pdf.savefig(fig)
                    plt.close(fig)

                self.gui.wykresy.plot_bar_chart()
                pdf.savefig(self.gui.wykresy.gui.bar_canvas.figure)

                self.gui.wykresy.plot_heatmap_visualization()
                pdf.savefig(self.gui.wykresy.gui.heatmap_canvas.figure)

                self.gui.wykresy.plot_motif_positions_visualization()
                pdf.savefig(self.gui.wykresy.gui.position_canvas.figure)

            QMessageBox.information(
                self.gui,
                "Eksport PDF",
                f"Raport zapisano do: {file_path}"
            )
            self.gui.log_output.append(f"Wyeksportowano raport PDF: {file_path}")

        except Exception as e:
            QMessageBox.critical(self.gui, "Błąd", f"Nie udało się zapisać PDF:\n{str(e)}")

    def summarize_statistics(self):
        if not hasattr(self.gui, "motif_df") or self.gui.motif_df.empty:
            return pd.DataFrame()

        pivot = self.gui.pivot_df.copy()
        sequences = [s["sequence"] for s in self.gui.sequences]

        summary = {
            "Liczba sekwencji": len(sequences),
            "Liczba aktywnych motywów": sum(m["active"] for m in self.gui.motifs),
            "Średnia liczba motywów na sekwencję": pivot.iloc[:, 1:].sum(axis=1).mean(),
            "Łączna liczba wszystkich motywów": pivot.iloc[:, 1:].sum().sum(),
            "Długość sekwencji min": min(len(seq) for seq in sequences),
            "Długość sekwencji max": max(len(seq) for seq in sequences),
            "Długość sekwencji średnia": sum(len(seq) for seq in sequences)/len(sequences)
        }

        df_summary = pd.DataFrame([summary])
        return df_summary

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = DNAAnalyzerGUI()
    window.show()
    sys.exit(app.exec())