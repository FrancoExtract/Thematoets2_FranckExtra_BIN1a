# Thematoets 2
# Franck Extra

try:
    import re
    import matplotlib.pyplot as plt
    import tkinter
    from tkinter import messagebox
except ModuleNotFoundError:
    print("Matplotlib module not found")  # It can't find the module for some reason...


def read_file():
    gff_file = open(gff, "r")
    # print(gff.readable())

    gff_file.close()


# def graph(x, y):
# plt.plot(x, y, 'g-')
# plt.title("Total counts of kinase")
# plt.xlabel("Chromosome")
# plt.ylabel("Kinase count")
# plt.show()

class GUI:
    def __init__(self):
        # Main window
        self.main_window = tkinter.Tk()

        # Top and bottom frames
        self.top_frame = tkinter.Frame(self.main_window)
        self.bottom_frame = tkinter.Frame(self.main_window)
        self.top_frame.pack()
        self.bottom_frame.pack()

        # Labels
        self.label1 = tkinter.Label(self.main_window, text="Hello top frame!")
        self.label1.pack()
        self.label2 = tkinter.Label(self.bottom_frame,
                                    text="Hello bottom frame!")
        self.label2.pack()

        # Buttons
        self.buttons = tkinter.Button(self.bottom_frame,  # Where the button is
                                      text="Click here lol",  # What text is shown
                                      command=self.do_something)  # What should happen when the button is pressed
        self.buttons.pack()
        self.quit_button = tkinter.Button(self.bottom_frame,
                                          text="Quit popup",
                                          command=self.main_window.destroy)
        self.quit_button.pack()

        # Shows main window
        tkinter.mainloop()

    def do_something(self):
        tkinter.messagebox.showinfo("Response",
                                    "Haha ur mom")


def main():
    gff = "GCF_000013425.1_ASM1342v1_genomic.gff"

    # GUI setup
    gui = GUI()

    # Prototype graph
    # graph(([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20]))


main()
