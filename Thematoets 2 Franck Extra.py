# Thematoets 2
# Franck Extra

import re
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox


class GFF:
    """
    This class defines the functions for reading and searching in
    the GFF file for the requested items within the file.

    Among these are:
    - The length of exons
    - The length of the gene
    - The protein ID associated with the gene
    - The accompanying CDS region
    """

    def __init__(self):
        """The init puts the starting position of all the functions to 0."""
        self.exon_amount = 0
        self.gene_length = 0
        self.protein_id = 0
        self.cds = 0

    def set_exons(self, exon):
        """This determines the amount of exons for the gene."""
        self.exon_amount = exon

    def get_exons(self):
        """This returns the determined amount of exons."""
        return self.exon_amount

    def set_genelength(self, start, stop):
        """
        This determines the length of the gene.
        This is done by subtracting the length of the stop with the start.

        Input:
        Whole gene = start --> str --> stop
        Calculation = stop - start
        """
        self.gene_length = int(stop) - int(start)

    def get_genelength(self):
        """This returns the calculated length of the gene."""
        return self.gene_length

    def set_proteinID(self, id):
        """This determines the ID for the protein."""
        self.protein_id = id

    def get_proteinID(self):
        """This returns the ID for the protein."""
        return self.protein_id

    def set_cds(self, code):
        """This determines and looks for the CDS region of the given gene."""
        self.cds = code

    def get_cds(self):
        """This returns the CDS region of the given gene."""
        return self.cds


def read_gff(gff):
    """
    This function reads the GFF file and adds the values both to
    their respective lists and objects within the GFF class.

    param gff: GFF file with exons, genes, and protein IDs
    """

    gff_file = open(gff, "r")  # Opens the GFF file
    if gff_file.readable():    # Boolean that prints if the file is readable
        print("The GFF file is readable and good to go.")
    else:
        print("The GFF file is unreadable.")
    print("")

    exon_list = []    # List for appending exons found in the file
    gene_length = []  # List for appending gene lengths found in the file

    # Reads the GFF file
    for line in gff_file:
        line = line.split("\t")  # Makes it a list
        # If NC_007795.1 in line
        if line[0] == "NC_007795.1":  # [0] is the first index in the column

            # If there's exon in the line
            if line[2] == "exon":  # 3rd column
                g1 = GFF()
                # Adds to set_exons
                g1.set_exons(int(line[4]) - int(line[3]))
                print(g1.get_exons())  # Prints the exon length
                exon_list.append(g1)   # Appends the exon list

            # If there's gene in the line
            if line[2] == "gene":
                g2 = GFF()
                # Adds to set_genelength
                g2.set_genelength(int(line[3]), int(line[4]))  # Gene length
                print(g2.get_genelength())  # Prints the length of the genes
                gene_length.append(g2)  # Appends to the gene length list

                # This block searches for protein IDs within the GFF file
                id_split = line[8].split(";")  # Splits the line into columns
                name = id_split[2]  # Defines the 2nd column as the ID name
                name.strip("Name=")  # Leaves only the ID name in the column
                g2.set_proteinID(name)  # Adds name to the protein ID object

    gff_file.close()


class GenBank_entry:
    """
    This class defines the objects in the GBFF file for:
    - Protein ID
    - Product (type of protein)
    - Protein sequence (translation)
    """

    def __init__(self, id, product, sequence):
        self.id = id
        self.product = product
        self.sequence = sequence


def read_gbff(gbff):
    """
    This function reads and retrieves the data for the objects
    defined in the GenBank class.

    :param gbff: GenBank file with the requested data
    """

    gbff_file = open(gbff, "r")  # Open the GBFF file
    if gbff_file.readable():  # Checks readability of the file
        print("The GBFF file is readable and good to go.")
    else:
        print("The GBFF file is unreadable.")
    print("")

    # Sets the default state of the values in the file lines to "False"
    gene = False
    CDS = False
    trans = False

    all_seq = []  # List for all the requested data in the file
    seq = []  # List for all available sequences in the file

    for line in gbff_file:
        # Calculates the length between the total line and empty spaces
        spaces = len(line) - len(line.lstrip(' '))

        if spaces == 5:  # Distance from "FEATURES" to "CDS, gene, locus" etc.
            if "gene" in line:  # Gene = True, when it's present in the line
                gene = True
                CDS = False
                trans = False
                seq = ''.join(seq)
                seq = seq.strip('/translation=')
                seq = seq.strip('"')
                if "id" in locals():
                    all_seq.append(GenBank_entry(id, product, seq))
                seq = []
            if "CDS" in line:  # CDS = True, when it's present in the line
                CDS = True
                gene = False

        if spaces == 21:  # Distance from "FEATURES" to the respective data
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


def serine_regex(sequence, regex_ser):
    """
    This function uses regex to find serine in the sequence.
    It returns "True" with a match; "False" is it doesn't.

    :param sequence: sequence for the serine kinase
    :param regex_ser: regex for serine
    :return: True/False
    """
    # Looks for a regex match in the sequence for serine
    match = re.search(regex_ser, sequence)
    regex_ser = "T.{2}[GC][NQ]SGS.[LIVM][FY]"
    if match:
        print("This is a serine sequence!")
        return True, regex_ser
    else:
        print("This is not a serine sequence.")
        return False, regex_ser


def histamine_regex(sequence, regex_his):
    """
    This function uses regex to find histamine in the sequence.
    It returns "True" with a match; "False" is it doesn't.

    :param sequence: sequence for the histamine kinase
    :param regex_his: regex for serine
    """
    # Looks for a regex match in the sequence for histamine
    match = re.search(regex_his, sequence)
    regex_his = "[ST]G[LIVMFYW]{3}[GN].{2}T[LIVM].T.{2}H"
    if match:
        print("This is a histamine sequence!")
        return True, regex_his
    else:
        print("This is not a histamine sequence.")
        return False, regex_his


def graph_prepare(x, y):
    """This function prepares a graph."""
    plt.plot(x, y, 'g-')
    plt.title("Total counts of kinase")
    plt.xlabel("Chromosome")
    plt.ylabel("Kinase count")
    plt.show()


class GUI:
    """This class sets up and gives output for the GUI. (Just shows results)"""
    def __init__(self):
        # Main window
        self.main_window = tk.Tk()             # Calls the window
        self.main_window.geometry("250x100")   # Dimensions of the window
        self.main_window.title('GFF Results')  # Title of the window

        # Top and bottom frames
        self.top_frame = tk.Frame(self.main_window)
        self.bottom_frame = tk.Frame(self.main_window)
        self.top_frame.pack()
        self.bottom_frame.pack()

        # Labels
        # self.label1 = tk.Label(self.main_window, text="""Click "Quit button" to terminate""")
        # self.label1.pack()
        self.label2 = tk.Label(self.bottom_frame, text="")
        self.label2.pack()

        """
        - The amount of exons
        - The length of the gene
        - The protein ID associated with the gene
        - The accompanying CDS region
        """

        # Buttons
        self.button1 = tk.Button(self.bottom_frame,  # Where the button is
                                 text="Click here for a list of the GFF items",  # What text is shown
                                 command=self.show_result)  # What should happen when the button is pressed
        self.button1.pack()

        self.quit_button = tk.Button(self.bottom_frame,
                                     text="Click here to quit the popup",
                                     command=self.main_window.destroy)
        self.quit_button.pack()

        # Shows main window
        tk.mainloop()

    @staticmethod
    def show_result():
        tk.messagebox.showinfo("GFF results per category",
                               "Check (one of) these boxes for the results")


def show_main_results():
    """
    This function globally shows the results retrieved from
    reading the GFF and GenBank file, and displays it per accession code.

    -
    """
    show_main_results()


def main():
    # GFF file
    gff = "C:/Users/franc/pythonProject/Weektaken/Course 2/Thematoets/GCF_000013425.1_ASM1342v1_genomic.gff.txt"
    read_gff(gff)

    # GenBank file
    gbff = "C:/Users/franc/pythonProject/Weektaken/Course 2/Thematoets/GCF_000013425.1_ASM1342v1_genomic.gbff.txt"
    read_gbff(gbff)
    
    # Serine regex
    serine_regex(sequence="", regex_ser="T.{2}[GC][NQ]SGS.[LIVM][FY]")
    
    # Histamine regex
    histamine_regex(sequence="", regex_his="[ST]G[LIVMFYW]{3}[GN].{2}T[LIVM].T.{2}H")
    
    # Prototype graph
    graph_prepare(([1, 2, 3, 4, 5]), ([2, 4, 6, 8, 10]))
    
    # GUI setup
    gui = GUI()
  
    return gff, gbff, gui


main()
