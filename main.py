"""
BigC
bigChem - tools for chemistry
6.14.2019
Python3
"""


import requests
import time
import string
from io import BytesIO
import sys
from PIL import ImageTk, Image
import curses




def getAll(formula_name):
    url_formula = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/property/MolecularFormula/TXT"
    url_structure = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/property/CanonicalSMILES/TXT"
    url_png = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/PNG"
    response = requests.get(url=url_png)
    data_png = response.content
    img = (Image.open(BytesIO(data_png)))
    response = requests.get(url=url_formula)
    data_formula = response.content
    response = requests.get(url=url_structure)
    data_structure = response.content
    return [data_formula, data_structure], img.show()


class createStructure:

    def __init__(self, formula, rname):
        self.formula = formula
        self.rname = rname
        formula = formula.decode("utf-8")
        formula = list(formula)

        while "\n" in formula:
            formula.remove("\n")

        self.formula = formula


    def molecular(self):
        rname = self.rname
        formula = self.formula
        structure = []
        returned_formula = formula
        returned_formula = "".join(returned_formula)

        for i in ["S", "C"]:
            if i in formula:
                x = formula.index(i)

                if formula[x+1] in [str(x) for x in range(101)]:

                    for _ in range(int(formula[x+1])):
                        y = list(i)
                        structure.extend(y)

        types = ["C", "H", "O", "N", "F", "S"]
        nums = {"S":0, "F":0, "C":0, "H": 0, "O" : 0, "N": 0}

        for i in formula:

            if i in types:
                z = formula.index(i)
                x = 1
                y = []

                if z == len(formula)-1:
                    break

                if len(formula[:z+x+1]) == len(formula):
                    for _ in formula[(z+1):]:
                        nums[i] += int(_)


                while formula[z+x] not in list(string.ascii_uppercase):


                    y.extend(formula[z+x:z+1+x])


                    x += 1

                    if z+x >= len(formula):
                        break


                if len(y) != 0:
                    if int("".join(y)) > 1:
                        nums[i] = int("".join(y))

                else:
                    nums[i] = 1

        count_in = 2
        count_out = len(structure) - 2

        for c, i in enumerate(structure):

            if int(nums["H"]) - 2 > -1 and int(nums["O"]) != 0 and count_out - 1 > -1:
                structure[c] = list("[" + str(i) + "-H-OH]")
                nums["H"] -= 2
                nums["O"] -= 1
                count_out -= 1


            elif int(nums["H"]) - 3 > -1 and int(nums["O"]) != 0 and count_in != 0:
                structure[c] = list("[" + str(i) + "H2OH]")
                nums["H"] -= 3
                nums["O"] -= 1
                count_in -= 1

            elif int(nums["H"]) != 0 and int(nums["O"]) != 0 and count_in != 0:
                structure[c] = list("[" + str(i) + "-H=O]")
                nums["H"] -= 1
                nums["O"] -= 1
                count_in -= 1

            elif int(nums["O"]) - 2 > -1 and count_in != 0:
                structure[c] = list("[" + str(i) + "=O=O]")
                nums["O"] -= 2
                count_in -= 1

            elif int(nums["H"]) - 3 > -1 and count_in != 0:
                structure[c] = list("[" + i + ("-H"*3) + "]")
                nums["H"] -= 3
                count_in -= 1

            elif nums["F"] - 5 > -1 and count_in != 0:
                structure[c] = list("[" + i + ("-F"*5) + "]")
                nums["F"] -= 5
                count_in -= 1

            elif nums["F"] - 3 > -1 and count_in != 0:
                structure[c] = list("[" + i + ("-F"*3) + "]")
                nums["F"] -= 3
                count_in -= 1

            elif nums["F"] - 2 > -1:
                structure[c] = list("[" + i + ("-F"*2) + "]")
                nums["F"] -= 2

            elif nums["H"] -2 > -1:
                structure[c] = list("[" + i + ("-H"*2) + "]")
                nums["H"] -= 2


        nums["C"] = 0
        nums["S"] = 0

        for i in list(nums.values()):
            if i != 0:
                if nums["O"] - 1 > -1 and nums["O"] - 1 > -1:
                    structure.append("[OH]")



        if len(structure) == 0:
            structure = formula

        struct = ""
        for i in structure:
            struct = struct + "".join(i)

        return struct



def checker(self, formula, rname):

    if "C" in formula:
        createStructure(formula, rname).molecular()


def compound():
    init = input(("\n"*100)+"Enter the name of the chemical structure you wish to use: ").lower()
    info = getAll(init)
    chem_formula = info[0]
    real_structure = chem_formula[1]
    chem_formula = chem_formula[0]

    ans = ((createStructure(chem_formula, init).molecular()))
    predicted_structure, real_structure = str(ans).replace("\n", ""), real_structure.decode("utf8")

    print("\nThe predicted structure is: {}.\nThe real structure is {}".format(predicted_structure, real_structure))
    print("The formula is: " + str(chem_formula.decode("utf-8").replace("\n", "")))
    print("The png structure is opened.")
    time.sleep(5)


menu = ["Periodic Table", 

    "[H ]                                                                  [He]",
    "[Li][Be]                                          [B ][C ][N ][O ][F ][Ne]",
    "[Na][Mg]                                          [Al][Si][P ][S ][Cl][Ar]",
    "[K ][Ca][Sc] [Ti][V ][Cr][Mn][Fe][Co][Ni][Cu][Zn] [Ga][Ge][As][Se][Br][Kr]",
    "[Rb][Sr][Y ] [Zr][Nb][Mo][Tc][Ru][Rh][Pd][Ag][Cd] [In][Sn][Sb][Te][I ][Xe]",
    "[Cs][Ba][La] [Hf][Ta][W ][Re][Os][Ir][Pt][Au][Hg] [Tl][Pb][Bi][Po][At][Rn]",
    "[Fr][Ra][Ac] [Rf][Db][Sg][Bh][Hs][Mt][Ds][Rg][Cn] [Nh][Fl][Mc][Lv][Ts][Og]",
    "                       ",
                 "[Ce][Pr][Nd][Pm][Sm][Eu][Gd][Tb][Dy][Ho][Er][Tm][Yb][Lu]",
                 "[Th][Pa][U ][Np][Pu][Am][Cm][Bk][Cf][Es][Fm][Md][No][Lr]",
    "               ",
        "Use the arrow keys to go over the periodic table! Press q to quit."
    ]


def printed(stdscr, current_idx):
    stdscr.clear()
    h, w = stdscr.getmaxyx()

    
    for idx, row in enumerate(menu):
        x = w//2 - len(row)//2
        y = h//2 - len(menu)//2 + idx
        if idx == current_idx:
            stdscr.attron(curses.color_pair(1))
            stdscr.addstr(y, x, row)
            stdscr.attroff(curses.color_pair(1))
        else: 
            stdscr.addstr(y, x, row)

    stdscr.refresh()



def periodic(stdscr):

    curses.curs_set(0)

    stdscr.keypad(True)
    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)

    current_row = 0
    current_pillar = 0
    printed(stdscr, current_row)

    while True:
        c = stdscr.getch()
        stdscr.clear()

        if c == ord("q"):
            break
        elif c == curses.KEY_UP and current_row > 0:
            current_row += 1
        elif c == curses.KEY_DOWN and current_row < len(menu)-1:
            current_row -= 1
        printed(stdscr, current_row)
        stdscr.refresh()

    curses.endwin()


def stoich():
    pass




#main loop
while True:
    print("\n"*75)
    print(
            """
                                          hN                                                      
   NMo         {}               ohysssyy  hN                                                       
   NMo`                       .mm:`   `.  hN_____                            
   NMdyooym.   ym   +ds+odm   yM/         hMysosdh    sho+oyn   .Nsyoodd+yoohm                     
   NMy    oM+  hM   Nh   +M   mM.         hN    NM   yM/    \:  /Ms`  -Md`  `mm                     
   NMo    :Ms  hM   hmyosds   yM/         hN    NM   dnsoooos   /M+   `My    mm                     
   NMd    hN   hM   Ny++/:    .dN+`   ./  hN    NM   os+        /M+   `My    mm                     
   Nh+syhhy    oh   dd+++odd   `/yhhhhyo  sh    Nh    yyhyyyy:  -h:   `ho    yy              
                   oM     hN                                                               
                    /ooooo+                                                     
            """
)
    print("                                  Compounds\n                               Periodic Table\n                         Stoichiometry (Coming soon!)")
    
    first = input("\nEnter what you want to use (type QUIT to exit): ")
    
    if first.lower() == "compounds" or first.lower() == "compound":
        compound()

    elif first.upper() == "QUIT":
        sys.exit()

    elif first.lower() == "periodic" or first.lower() == "periodic table":
        curses.wrapper(periodic)

