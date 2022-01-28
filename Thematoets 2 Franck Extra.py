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

    def set_exon_length(self, exon):
        """This determines the amount of exons for the gene."""
        self.exon_amount = exon

    def get_exon_length(self):
        """This returns the determined amount of exons."""
        return self.exon_amount

    def set_gene_length(self, start, stop):
        """
        This determines the length of the gene.
        This is done by subtracting the length of the stop with the start.

        Input:
        Whole gene = start --> str --> stop
        Calculation = stop - start
        """
        self.gene_length = int(stop) - int(start)

    def get_gene_length(self):
        """This returns the calculated length of the gene."""
        return self.gene_length

    def set_protein_ID(self, id):
        """This determines the ID for the protein."""
        self.protein_id = id

    def get_protein_ID(self):
        """This returns the ID for the protein."""
        return self.protein_id

    def set_CDS(self, code):
        """This determines and looks for the CDS region of the given gene."""
        self.cds = code

    def get_CDS(self):
        """This returns the CDS region of the given gene."""
        return self.cds


def read_gff(gff):
    """
    This function reads the GFF file and adds the values both to
    their respective lists and objects within the GFF class.

    Among these are:
    - The length of exons
    - The length of genes
    - The identification tag of proteins
    - The CDS region

    param gff: exon, gene, protein ID
    """

    global gff_file  # Assignment for "gff_file" in line 100

    try:
        gff_file = open(gff, "r")  # Opens the GFF file
        if gff_file.readable():    # Boolean that prints if the file is readable
            print("The GFF file is readable and good to go.")
        else:
            print("The GFF file is unreadable.")
        print("")
    except FileNotFoundError as err:  # If the file is not in the same folder
        print(err)
        print("Make sure the file is in the same folder, and try again.")

    exon_list = []    # List for exons found in the file
    gene_length = []  # List for gene lengths found in the file

    # Reads the GFF file
    for line in gff_file:
        line = line.split("\t")  # Makes it a list
        # If NC_007795.1 in line
        if line[0] == "NC_007795.1":  # [0] is the first index in the column

            # If there's "exon" in the third column of the line
            if line[2] == "exon":  # 3rd column
                g1 = GFF()  # g1 is designated for retrieving exons
                g1.set_exon_length(
                    int(line[4]) - int(line[3]))  # Adds to set_exons
                print(g1.get_exon_length())  # Prints the exon length
                exon_list.append(g1)  # Appends the exon list

            # If there's "gene" in the third column of the line
            if line[2] == "gene":
                g2 = GFF()  # g2 is designated for genes and protein IDs
                g2.set_gene_length(int(line[3]), int(line[4]))  # Gene length
                gene_length.append(g2)  # Appends to the gene length list

                # This block searches for protein IDs within the GFF file
                id_split = line[8].split(";")  # Splits the line into columns
                name = id_split[2]  # Defines the 2nd column as the ID name
                name.strip("Name=")  # Leaves only the ID name in the column
                g2.set_protein_ID(name)  # Adds name to the protein ID object

    gff_file.close()  # Closes the GFF file after reading
    return exon_list, gene_length, GFF


class GenBank_entries:
    """
    This class defines the objects in the GBFF file for:
    - Protein ID
    - Product (type of protein)
    - Protein sequence (translation)
    """

    def __init__(self, id, product, sequence):
        self.id = id  # Protein ID
        self.product = product  # Name of the protein
        self.sequence = sequence  # Protein sequence


def read_gbff(gbff):
    """
    This function reads and retrieves the data for the objects
    defined in the GenBank class.

    :param gbff: GenBank file with the requested data
    """

    global id, product  # Assignment for "id" and "product" in line 180

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
                seq = ''.join(seq)  # Joins any spaces together in the line
                seq = seq.strip('/translation=')  # Removes transl. from line
                seq = seq.strip('"')  # Strips out any space in the line
                if "id" in locals():  # Adds ID to entry if present in line
                    all_seq.append(GenBank_entries(id, product, seq))
                seq = []
            if "CDS" in line:  # CDS = True, when it's present in the line
                CDS = True
                gene = False

        if spaces == 21:  # Distance from "FEATURES" to the respective data
            if CDS is True:
                if "/protein_id" in line:
                    id = line.replace('/protein_id=', '')  # Replace with space
                    id = id.replace('"', '')  # Replace 2 spaces with 1 space
                    id = id.strip(' ')  # Removes the 2 combined spaces
                if '/product' in line:
                    p = line.replace('/product=', '')  # Product becomes space
                    p = p.replace('"', '')
                    p = p.strip(' ')
                    product = p
                if '/translation' in line:
                    trans = True
                if trans and not gene:
                    line = line.strip('\n')  # Strips out any newlines
                    line = line.strip(' ')  # Strips out any 2 spaces
                    seq.append(line)  # Appends line to the sequence list

    for x in all_seq:  # "x" returns any value within all_seq
        print(x.id)

    gbff_file.close()  # Closes the GBFF file after reading

    # "all_seq" reads the whole GenBank file
    return all_seq, seq, id, trans, gene, CDS, product


def serine_regex(sequence, regex_ser, all_seq=None):
    """
    This function uses regular expression to find serine in the sequence.
    It returns "True" if it matches; "False" if it doesn't.

    :param sequence: sequence for the serine kinase
    :param regex_ser: regular expression for serine
    :param all_seq: extra parameter for the GenBank information
    :return: True/False
    """

    # Looks for a regex match in the sequence for serine
    match = re.search(regex_ser, sequence)  # Compares regex with the sequence
    regex_ser = "T.{2}[GC][NQ]SGS.[LIVM][FY]"  # Serine regex
    for x in all_seq:
        sequence = x.seq
    if match:
        print("This is a serine sequence!")
        return True, regex_ser, sequence
    else:
        print("This is not a serine sequence.")
        return False, regex_ser, sequence


def histamine_regex(sequence, regex_his, all_seq=None):
    """
    This function uses regular expression to find histamine in the sequence.
    It returns "True" if it matches; "False" if it doesn't.

    :param sequence: sequence for the histamine kinase
    :param regex_his: regular expression for histamine
    :param all_seq: extra parameter for the GenBank information
    """

    # Looks for a regex match in the sequence for histamine
    match = re.search(regex_his, sequence)  # Compares regex with the sequence
    regex_his = "[ST]G[LIVMFYW]{3}[GN].{2}T[LIVM].T.{2}H"  # Histamine regex
    for x in all_seq:
        sequence = x.seq
    if match:
        print("This is a histamine sequence!")
        return True, regex_his, sequence
    else:
        print("This is not a histamine sequence.")
        return False, regex_his, sequence


def graph_prepare(x, y):
    """
    This function prepares a graph. (I don't yet know what I'm using it for)

    :param x: chromosome data on the x-axis
    :param y: exon data on the y-axis
    """

    plt.plot(x, y, 'g-')
    plt.title("Total counts of exons per chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Exon count")
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
        self.label1 = tk.Label(self.bottom_frame, text="")  # Top line is blank
        self.label1.pack()

        # Buttons
        self.button1 = tk.Button(self.bottom_frame,  # Button location
                                 text="Click here for a list of the GFF items",  # Text shown
                                 command=self.show_result)  # Action when button is pressed
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
                               "The results are shown in the terminal :)")


def show_main_results(read_gff, read_gbff):
    """
    This function globally shows the results retrieved from
    reading the GFF and GenBank file, and displays it per accession code.

    :param read_gff: items pertaining to the GFF file
    :param read_gbff: items pertaining to the GBFF file
    """

    print(read_gff)
    print(read_gbff)
    print(GFF)


def main():
    # GFF file
    gff = "GCF_000013425.1_ASM1342v1_genomic.gff.txt"
    read_gff(gff)

    # GenBank file
    gbff = "GCF_000013425.1_ASM1342v1_genomic.gbff.txt"
    read_gbff(gbff)
    
    # Serine regex
    x = all_seq.trans
    # serine_regex(sequence=x, regex_ser="T.{2}[GC][NQ]SGS.[LIVM][FY]")

    # Histamine regex
    x = all_seq.trans
    # histamine_regex(sequence=x, regex_his="[ST]G[LIVMFYW]{3}[GN].{2}T[LIVM].T.{2}H")
    
    # Prototype graph
    # graph_prepare(([1, 2, 3, 4, 5]), ([2, 4, 6, 8, 10]))

    # GUI setup
    # gui = GUI()

    # Shows the main results
    # show_main_results(read_gff, read_gbff)

    return gff, gbff, x,  # gui


main()
