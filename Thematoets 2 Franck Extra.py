# Thematoets 2
# Franck Extra

import regex as re
import matplotlib.pyplot as plt
import tkinter as tk


class GFF3:
    """
    Deze class definieert de functies voor het inlezen van de GFF3 file.

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

    def set_exonen(self, ex):
        """
        Dit bepaalt het aantal exonen voor het gen.
        """
        self.aantal_exonen = ex

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

    def set_accessiecode(self, accessiecode):
        """
        Dit bepaalt en zoekt naar de accessiecode van het opgegeven gen.
        """
        self.accessiecode = accessiecode

    def get_accessiecode(self):
        """
        Dit returnt de accessiecode van het opgegeven gen
        """
        return self.accessiecode


if __name__ == "__main__":
    # De verschillende variabelen
    regex_ser = "T-x(2)-[GC]-[NQ]-S-G-S-x-[LIVM]-[FY]"  # Serine regex
    regex_his = "[ST]-G-[LIVMFYW](3)-[GN]-x(2)-T-[LIVM]-x-T-x(2)-H"  # Histamine regex
    gff3 = "GCF_000013425.1_ASM1342v1_genomic.gff"
