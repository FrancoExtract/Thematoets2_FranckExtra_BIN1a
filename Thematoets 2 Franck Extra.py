# Thematoets 2
# Franck Extra

import re
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox


class GFF:
    """
    Deze class definieert de functies voor het inlezen van
    en zoeken in de GFF file naar gevraagde items binnen de file.

    Dit zijn onder andere:
    - Het aantal exonen
    - De lengte van het gen
    - De benodigde accessiecode van het gen
    """

    def __init__(self):
        # De init zet de beginstand van de functies allemaal naar 0
        self.aantal_exonen = 0
        self.lengte_gen = 0
        self.chromosoom = 0
        self.accessiecode = 0

    def set_exonen(self, exon):
        """Dit bepaalt het aantal exonen voor het gen."""
        self.aantal_exonen = exon

    def get_exonen(self):
        """Dit returnt het bepaalde aantal exonen."""
        return self.aantal_exonen

    def set_lengtegen(self, start, stop):
        """
        Dit bepaalt de lengte van het gen.
        Dit gebeurt door de lengte van de stop min start te berekenen.

        Input:
        Gehele gen = start - str - stop
        Berekening = stop - start
        """
        self.lengte_gen = int(stop) - int(start)

    def get_lengtegen(self):
        """Dit returnt de berekende lengte van het gen."""
        return self.lengte_gen

    def set_proteinID(self, chr):
        """Dit bepaalt het ID voor het eiwit."""
        self.chromosoom = chr

    def get_proteinID(self):
        """Dit returnt het ID voor het eiwit."""
        return self.chromosoom

    def set_cds(self, code):
        """Dit bepaalt en zoekt naar de CDS region van het opgegeven gen."""
        self.accessiecode = code

    def get_cds(self):
        """Dit returnt de CDS region van het opgegeven gen."""
        return self.accessiecode


def read_gff(gff):
    gff_file = open(gff, "r")
    # print(gff_file.readable())

    exon_list = []
    gene_length = []

    for line in gff_file:
        line = line.split("\t")  # Makes it a list
        # If NC_007795.1 in line
        if line[0] == "NC_007795.1":  # [0] is the first index in the column
            # If line isn't empty
            if line[2] == "exon":  # 3rd column
                g1 = GFF()
                # Adds to set_exons
                g1.set_exonen(int(line[4]) - int(line[3]))
                print(g1.get_exonen())  # Prints the length of the exons
                exon_list.append(g1)  # Adds the g1 object to the g1 list

            if line[2] == "gene":
                g2 = GFF()
                # Adds to set_genelength
                g2.set_lengtegen(int(line[3]), int(line[
                                                       4]))  # Parameters for the calculation were already set
                print(g2.get_lengtegen())
                gene_length.append(g2)

                id_split = line[8].split(";")
                name = id_split[2]
                name.strip("Name=")
                g2.set_proteinID(name)

    gff_file.close()


class genebank_entry:
    def __init__(self, id, product, sequence):
        self.id = id
        self.product = product
        self.sequence = sequence


def read_gbff(gbff):
    gbff_file = open(gbff, "r")
    # print(gbff_file.readable())

    gene = False
    CDS = False
    trans = False

    all_seq = []
    seq = []

    for line in gbff_file:
        spaces = len(line) - len(line.lstrip(' '))

        if spaces == 5:
            if "gene" in line:
                gene = True
                CDS = False
                trans = False
                seq = ''.join(seq)
                seq = seq.strip('/translation=')
                seq = seq.strip('"')
                if "id" in locals():
                    all_seq.append(genebank_entry(id, product, seq))
                seq = []
            if "CDS" in line:
                CDS = True
                gene = False

        if spaces == 21:
            if CDS is True:
                if "/protein_id" in line:
                    id = line.replace('/protein_id=', '')
                    id = id.replace('"', '')
                    id = id.strip(' ')

                if '/product' in line:
                    p = line.replace('/product=', '')
                    p = p.replace('"', '')
                    p = p.strip(' ')
                    product = p
                if '/translation' in line:
                    trans = True
                if trans and not gene:
                    line = line.strip('\n')
                    line = line.strip(' ')
                    seq.append(line)

    for x in all_seq:
        print(x.id)

    gbff_file.close()

    # all_seq reads the whole GenBank file
    return all_seq


def set_serine(sequence, regex_ser):
    """
    Deze functie gebruikt regex om serine in de sequentie te vinden.

    :param sequence: sequentie voor de kinase
    :param regex_ser: regex voor serine
    :return: True/False
    """
    # Zoekt naar een regex-match in de sequentie voor serine
    regex_ser = 'T.{2}[GC][NQ]SGS.[LIVM][FY]'
    match = re.search(regex_ser, sequence)
    if match:
        return True
    else:
        return False


def histamine_regex(sequence, regex_his):
    """
    Deze functie gebruikt regex om histamine in de sequentie te vinden.
    Als het matcht, returnt het "True"; zo niet, dan returnt het "False".

    :param sequence: sequentie voor de kinase
    :param regex_his: regex voor serine
    """
    # Zoekt naar een regex-match in de sequentie voor histamine
    regex_his = "[ST]G[LIVMFYW]{3}[GN].{2}T[LIVM].T.{2}H"
    match = re.search(regex_his, sequence)
    if match:
        return True
    else:
        return False


# def graph(x, y):
# plt.plot(x, y, 'g-')
# plt.title("Total counts of kinase")
# plt.xlabel("Chromosome")
# plt.ylabel("Kinase count")
# plt.show()


class GUI:
    def __init__(self):
        # Main window
        self.main_window = tk.Tk()

        # Top and bottom frames
        self.top_frame = tk.Frame(self.main_window)
        self.bottom_frame = tk.Frame(self.main_window)
        self.top_frame.pack()
        self.bottom_frame.pack()

        # Labels
        self.label1 = tk.Label(self.main_window, text="Mfw")
        self.label1.pack()
        self.label2 = tk.Label(self.bottom_frame, text="Bottom text")
        self.label2.pack()

        # Buttons
        self.buttons = tk.Button(self.bottom_frame,  # Where the button is
                                 text="Click here lol",  # What text is shown
                                 command=self.do_something)  # What should happen when the button is pressed
        self.buttons.pack()
        self.quit_button = tk.Button(self.bottom_frame,
                                     text="Quit popup",
                                     command=self.main_window.destroy)
        self.quit_button.pack()

        # Shows main window
        tk.mainloop()

    def do_something(self):
        tk.messagebox.showinfo("Response",
                               "Haha ur mom")


def main():
    # GFF file
    gff = "C:/Users/franc/pythonProject/Weektaken/Course 2/Thematoets/GCF_000013425.1_ASM1342v1_genomic.gff.txt"
    read_gff(gff)

    # GenBank file
    gbff = "C:/Users/franc/pythonProject/Weektaken/Course 2/Thematoets/GCF_000013425.1_ASM1342v1_genomic.gbff.txt"
    read_gbff(gbff)

    # GUI setup
    gui = GUI()

    # Prototype graph
    # graph(([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20]))

    return gff


main()
