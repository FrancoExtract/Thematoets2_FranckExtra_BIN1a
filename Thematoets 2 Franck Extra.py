# Thematoets 2
# Franck Extra

# Import the needed modules for this Python file
try:
    import re
    import matplotlib.pyplot as plt
    import tkinter as tk
    from tkinter import messagebox
except ModuleNotFoundError as err1:  # If a certain module is not installed
    print(err1)


class GFF:
    """
    This class defines the functions for reading and searching in
    the GFF file for the requested items within the file.

    Among these are:
    - The length of exons
    - The length of the gene
    - The protein ID associated with the gene
    - The accompanying CDS length
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

    def set_protein_id(self, id):
        """This determines the ID for the protein."""
        self.protein_id = id

    def get_protein_id(self):
        """This returns the ID for the protein."""
        return self.protein_id

    def set_cds_length(self, code):
        """This determines and looks for the CDS length of the given gene."""
        self.cds = code

    def get_cds_length(self):
        """This returns the CDS length of the given gene."""
        return self.cds


def read_gff(gff):
    """
    This function reads the GFF file and adds the values both to
    their respective lists and objects within the GFF class.

    Among these are:
    - The length of exons
    - The length of genes
    - The identification tag of proteins
    - The CDS length

    :param gff: exon, gene, protein ID from GFF file
    :return: exon_list, gene_length, CDS
    """

    print("GFF Data:")
    print("")
    try:
        gff_file = open(gff, "r")  # Opens the GFF file
        if gff_file.readable():  # Boolean that checks file readability
            print("(The GFF file is readable and good to go)")
        print("")
    except FileNotFoundError as err2:  # If the file is not in the same folder
        print(err2)
        print("Make sure the file is in the same folder, and try again.")
        print("")

    exon_list = []  # List for exon lengths found in the file
    gene_length = []  # List for gene lengths found in the file
    CDS_length = []  # List for the CDS lengths found in the file

    # Reads the GFF file
    for line in gff_file:
        line = line.split("\t")  # Makes it a list

        # If "NC_007795.1" in line
        if line[0] == "NC_007795.1":  # [0] is the first index in the column

            # If there's "exon" in the third column of the line
            if line[2] == "exon":  # 3rd column
                g1 = GFF()  # g1 is designated for retrieving exons
                g1.set_exon_length(
                    int(line[4]) - int(line[3]))  # Adds to set_exons
                exon_list.append(g1)  # Appends the exon list

            # If there's "gene" in the third column of the line
            if line[2] == "gene":
                g2 = GFF()  # g2 is designated for genes and protein IDs
                g2.set_gene_length(int(line[3]), int(line[4]))  # Gene length
                gene_length.append(g2)  # Appends to the gene length list

                # This block searches for protein IDs within the GFF file
                id_split = line[8].split(";")  # Splits the line into columns
                name = id_split[2]  # Defines the 2nd column as the ID name
                name = name.replace("Name=", "")  # Leaves only the ID name
                g2.set_protein_id(name)  # Adds name to the protein ID object

            # If there's "CDS" in the third column of the line
            if line[2] == "CDS":
                g3 = GFF()  # g3 is designated for the CDS length
                g3.set_cds_length(
                    int(line[4]) - int(line[3]))  # Adds to set_CDS
                CDS_length.append(g3)  # Appends to the CDS list

    gff_file.close()  # Closes the GFF file after reading
    return exon_list, gene_length, CDS_length


class GenBank_entries:
    """
    This class defines the objects (entries) in the GBFF file for:
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
    :return: all_data
    """

    print("GBFF Data:")
    print("")

    try:
        gbff_file = open(gbff, "r")  # Open the GBFF file
        if gbff_file.readable():  # Boolean that checks file readability
            print("(The GBFF file is readable and good to go)")
        print("")
    except FileNotFoundError as err3:  # If the file is not in the same folder
        print(err3)
        print("Make sure the file is in the same folder, and try again.")
    print("")

    # Sets the default state of the values in the file lines to "False"
    gene = False
    CDS = False
    trans = False

    all_data = []  # List for all the requested data in the file
    seq = []  # List for all available sequences (trans) in the file

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
                    all_data.append(GenBank_entries(id, product, seq))
                seq = []
            if "CDS" in line:  # CDS = True, when it's present in the line
                CDS = True
                gene = False

        if spaces == 21:  # Distance from "FEATURES" to the respective data
            if CDS is True:
                if "/protein_id" in line:
                    id = line.replace('/protein_id=', '')  # Leaves blank space
                    id = id.replace('"', '')  # Replaces quote with blank space
                    id = id.strip(' ')  # Strips out any leftover spaces
                if '/product' in line:
                    p = line.replace('/product=', '')  # Leaves blank space
                    p = p.replace('"', '')  # Replace quote with blank space
                    p = p.strip(' ')  # Strips out any leftover spaces
                    product = p
                if '/translation' in line:
                    trans = True
                if trans and not gene:
                    line = line.strip('\n')  # Strips out any newlines
                    line = line.strip(' ')  # Strips out any leftover space
                    seq.append(line)  # Appends line to the sequence list

    gbff_file.close()  # Closes the GBFF file after reading

    return all_data  # "all_data" reads the whole GenBank file


def serine_regex(sequence, regex_ser):
    """
    This function uses regular expression to find serine in the sequence.
    It returns "True" if it matches; "False" if it doesn't.

    :param sequence: sequence for the serine kinase
    :param regex_ser: regular expression for serine
    :return: True/False
    """

    # Looks for a regex match in the sequence for serine
    match = re.search(sequence, regex_ser)  # Compares regex with the sequence
    if match:
        return True
    else:
        return False


def histamine_regex(sequence, regex_his):
    """
    This function uses regular expression to find histamine in the sequence.
    It returns "True" if it matches; "False" if it doesn't.

    :param sequence: sequence for the histamine kinase
    :param regex_his: regular expression for histamine
    """

    # Looks for a regex match in the sequence for histamine
    match = re.search(sequence, regex_his)  # Compares regex with the sequence
    if match:
        return True
    else:
        return False


def graph_prepare(x, y):
    """
    This function prepares a graph.
    (I don't yet know how I will implement the use for it though)

    :param x: chromosome data on the x-axis
    :param y: exon data on the y-axis
    """

    plt.plot(x, y, 'g-')  # Plots out the axes and makes the line green
    plt.title("Total counts of exons per chromosome")  # Title for the graph
    plt.xlabel("Chromosome")  # Title for the x-axis
    plt.ylabel("Exon count")  # Title for the y-axis
    plt.show()  # Makes the graph show up when running the program


class GUI:
    """This class sets up and gives output for the GUI. (Just shows results)"""

    def __init__(self):
        # Main window
        self.main_window = tk.Tk()  # Calls the window
        self.main_window.geometry("300x100")  # Dimensions of the window
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
                                 text="Click here to see where the results "
                                      "will be printed",  # Text shown
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
        tk.messagebox.showinfo("GFF & GenBank Results",
                               "The results will show in the terminal "
                               "after closing the popup :)")


def main():
    """
    Here, the main functions are called and most of the variables within this
    project declared.

    :return: gff, gbff, gui
    """

    # GUI setup
    gui = GUI()

    # Prototype graph
    graph_prepare(([1, 2, 3, 4, 5]), ([2, 4, 6, 8, 10]))

    # GFF file
    gff = "GCF_000013425.1_ASM1342v1_genomic.gff.txt"
    GFF_list = read_gff(gff)  # GFF returns everything in "read_gff"

    # GFF prints
    for nested_list in GFF_list:  # List in the GFF list
        for obj in nested_list:  # Object within the list's list
            print("Exon length: ", obj.get_exon_length(), "|",
                  "Gene length: ", obj.get_gene_length(), "|",
                  "Protein ID: ", obj.get_protein_id(), "|",
                  "CDS length: ", obj.get_cds_length())
    print("-" * 80)
    print("")

    # GenBank file
    gbff = "GCF_000013425.1_ASM1342v1_genomic.gbff.txt"
    all_data = read_gbff(gbff)

    # Serine regex prints
    print("These are serine sequences:")
    for x in all_data:
        ser_match = serine_regex(sequence=x.sequence,
                                 regex_ser="T.{2}[GC][NQ]SGS.[LIVM][FY]")
        print(ser_match)
        if ser_match:
            print(x.sequence)
    print("")

    # Histamine regex prints
    print("These are histamine sequences:")
    for x in all_data:
        # Histamine regex
        his_match = histamine_regex(sequence=x.sequence,
                                    regex_his="[ST]G[LIVMFYW]{3}[GN].{2}T[LIVM].T.{2}H")
        print(his_match)
        if his_match:
            print(x.sequence)
    print("")
    print("""For some reason, regex *does* match, but it prints only "False",
    and "True"; I don't know why...""")
    print("-" * 80)
    print("")

    return gff, gbff, gui


main()
