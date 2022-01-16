# Thematoets 2
# Franck Extra

import regex as re
import matplotlib.pyplot as plt
import tkinter as tk


class GFF:
    """
    Deze class definieert de functies voor het inlezen van
    en zoeken in de GFF file.

    Dit zijn onder andere:
    - Het aantal exonen
    - De lengte van het gen
    - Het chromosoom voor het gen
    - De benodigde accessiecode van het gen
    """

    def __init__(self):
        # De init zet de beginstand van de functies allemaal naar 0
        self.aantal_exonen = 0
        self.lengte_gen = 0
        self.chromosoom = 0
        self.accessiecode = 0

    def set_exonen(self, exon):
        """
        Dit bepaalt het aantal exonen voor het gen.
        """
        self.aantal_exonen = exon

    def get_exonen(self):
        """
        Dit returnt het bepaalde aantal exonen.
        """
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
        """
        Dit returnt de berekende lengte van het gen.
        """
        return self.lengte_gen

    def set_chromosoom(self, chr):
        """
        Dit bepaalt het chromosoom voor het opgegeven gen.
        """
        self.chromosoom = chr

    def get_chromosoom(self):
        """
        Dit returnt het chromosoom van het opgegeven gen.
        """
        return self.chromosoom

    def set_accessiecode(self, code):
        """
        Dit bepaalt en zoekt naar de accessiecode van het opgegeven gen.
        """
        self.accessiecode = code

    def get_accessiecode(self):
        """
        Dit returnt de accessiecode van het opgegeven gen.
        """
        return self.accessiecode


def serine_regex(sequence, regex_ser):
    """
    Deze functie gebruikt regex om serine in de sequentie te vinden.
    Als het matcht, returnt het "True"; zo niet, dan returnt het "False".

    :param sequence: sequentie voor de kinase
    :param regex_ser: regex voor serine
    :return: True/False
    """
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
    :return: True/False
    """
    match = re.search(regex_his, sequence)
    if match:
        return True
    else:
        return False


def gff_list(gff, fasta_dict):
    """
    Dit maakt een lijst met daarin alle FASTA-items,
    plus informatie uit de GFF-class:
    - Exonen
    - Genlengte
    - Chromosoom
    - Accessiecode van het gen

    :param gff: str - GFF file
    :param fasta_dict: dict - fasta dictionary {header: seq}
    :return: lijst met GFF-items
    """
    gff_entries = []
    entry = ""
    # Opent en leest het bestand in
    try:
        with open(gff) as inFile:
            for line in inFile:
                # Als er "mRNA" in de regel te vinden is
                if re.search("mRNA", line):
                    # Als de entry níét leeg is
                    if entry != "":
                        entry.set_exonen(exonen)
                        if entry.get_accessiecode() in fasta_dict.keys():
                            gff_entries.append(entry)
                            exonen = 0

                        # Dit maakt een nieuwe entry aan
                        entry = GFF()
                        entry.set_chromosoom(line.split()[0].split("Chr")[1])
                        entry.set_lengtegen(line.split()[3], line.split()[4])
                        entry.set_accessiecode(line.split()[8].split(";")[0] \
                                               .split("=")[1])

                        # Dit haalt het chromosoom, start-stop en
                        # de accessiecode uit de regel
                        entry.set_chromosoom(line.split()[0].split("Chr")[1])
                        entry.set_lengtegen(line.split()[3], line.split()[4])
                        entry.set_accessiecode(line.split()[8].split(";")[0] \
                                               .split("=")[1])

                # Telt +1 met het aantal exonen dat wordt gevonden
                elif re.search("exon", line).lower():
                    exonen += 1

    except FileNotFoundError:  # Als het bestand niet in dezelfde map staat
        print("De gevraagde file is niet aanwezig:", "GCF_000013425.1_ASM1342v1_genomic.gff")
        gff = input("Geef een nieuw bestand: ")
        gff_list(gff, fasta_dict)

        return gff_entries


if __name__ == "__main__":
    # De verschillende variabelen
    regex_ser = "T-x(2)-[GC]-[NQ]-S-G-S-x-[LIVM]-[FY]"  # Serine regex
    regex_his = "[ST]-G-[LIVMFYW](3)-[GN]-x(2)-T-[LIVM]-x-T-x(2)-H"  # Histamine regex
    gff = "GCF_000013425.1_ASM1342v1_genomic.gff"
    gbff_file = "GCF_000013425.1_ASM1342v1_genomic.gbff"
