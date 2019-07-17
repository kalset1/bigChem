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
import math
from mendeleev import element, IonicRadius
import string

menu = [
    "[H ]", "[He]","[Li]", "[Be]", "[B ]", "[C ]", "[N ]", "[O ]", "[F ]", "[Ne]",
    "[Na]", "[Mg]", "[Al]", "[Si]", "[P ]", "[S ]", "[Cl]", "[Ar]",
    "[K ]", "[Ca]", "[Sc]", "[Ti]", "[V ]", "[Cr]", "[Mn]", "[Fe]", "[Co]", "[Ni]", "[Cu]", "[Zn]", "[Ga]", "[Ge]", "[As]", "[Se]", "[Br]", "[Kr]",
    "[Rb]", "[Sr]", "[Y ]", "[Zr]", "[Nb]", "[Mo]", "[Tc]", "[Ru]", "[Rh]", "[Pd]", "[Ag]", "[Cd]", "[In]", "[Sn]", "[Sb]", "[Te]", "[I ]", "[Xe]",
    "[Cs]", "[Ba]", "[La]", "[Hf]", "[Ta]", "[W ]", "[Re]", "[Os]", "[Ir]", "[Pt]", "[Au]", "[Hg]", "[Tl]", "[Pb]", "[Bi]", "[Po]", "[At]", "[Rn]",
    "[Fr]", "[Ra]", "[Ac]", "[Rf]", "[Db]", "[Sg]", "[Bh]", "[Hs]", "[Mt]", "[Ds]", "[Rg]", "[Cn]", "[Nh]", "[Fl]", "[Mc]", "[Lv]", "[Ts]", "[Og]",
    "[Ce]", "[Pr]", "[Nd]", "[Pm]", "[Sm]", "[Eu]", "[Gd]", "[Tb]", "[Dy]", "[Ho]", "[Er]", "[Tm]", "[Yb]", "[Lu]",
    "[Th]", "[Pa]", "[U ]", "[Np]", "[Pu]", "[Am]", "[Cm]", "[Bk]", "[Cf]", "[Es]", "[Fm]", "[Md]", "[No]", "[Lr]"
    ]

elements = ["Hydrogen",'Helium','Lithium','Beryllium','Boron','Carbon','Nitrogen','Oxygen','Fluorine','Neon','Sodium','Magnesium','Aluminum','Silicon','Phosphorus','Sulfur','Chlorine',
'Argon','Potassium','Calcium','Scandium','Titanium','Vanadium','Chromium','Manganese','Iron','Cobalt','Nickel','Copper','Zinc','Gallium','Germanium','Arsenic',
'Selenium','Bromine','Krypton','Rubidium','Strontium','Yttrium','Zirconium','Niobium','Molybdenum','Technetium','Ruthenium','Rhodium','Palladium',
'Silver','Cadmium','Indium','Tin','Antimony','Tellurium','Iodine','Xenon','Cesium','Barium',
'Lanthanum','Cerium','Praseodymium','Neodymium','Promethium','Samarium','Europium','Gadolinium','Terbium','Dysprosium','Holmium','Erbium','Thulium','Ytterbium','Lutetium','Hafnium','Tantalum','Tungsten','Rhenium','Osmium','Iridium','Platinum','Gold','Mercury','Thallium','Lead','Bismuth','Polonium','Astatine','Radon','Francium','Radium','Actinium','Thorium','Protactinium','Uranium','Neptunium','Plutonium','Americium',
'Curium','Berkelium','Californium','Einsteinium','Fermium','Mendelevium','Nobelium','Lawrencium','Rutherfordium','Dubnium','Seaborgium','Bohrium',
'Hassium','Meitnerium','Darmstadtium','Roentgenium','Copernicium',"Nihonium","Flerovium","Moscovium","Livermorium","Tennessine","Oganesson"]

values = []

for i in elements:
    x = element(i)
    symbol = x.symbol
    atom_number = x.atomic_number
    atom_weight = x.mass
    atom_radius = x.atomic_radius
    config = x.econf
    y = [symbol, atom_number, atom_weight, config]
    values.append(y)

element_dict = {}

for i in range(len(elements)):
    element_dict[elements[i]] = values[i]

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


    def ionic(self):
        charges = []

        formula = self.formula
        rname = self.rname
        structure = "[" + "".join(formula) + "]"

        formula[0], formula[1] = (str(formula[0]) + str(formula[1])), (str(formula[2]) + str(formula[3]))

        for i in formula:
            if i not in str(range(101)):
                z = element(i)
                z = z.ionic_radii
                z = list(z)
                
                for i in z:
                    if "charge=" in str(i):
                        i = str(i).replace("charge=", "")
                        i = list(i)
                        i = i[:i.index(",")]
                        i = "".join(i)
                        i = str(i).strip()
                        charges.append(int(i))

        return [structure, charges]
        

    def molecular(self):
        rname = self.rname
        formula = self.formula
        structure = []
        returned_formula = formula
        returned_formula = "".join(returned_formula)

        if "l" in formula and "C" in formula:
            if formula.index("C") - formula.index("l") == -1:
                cindex = formula.index("C")
                del formula[cindex] 
                del formula[formula.index("l")]
                formula[cindex] = "Cl"

        for i in ["S", "C"]:
            if i in formula:
                x = formula.index(i)

                if formula[x+1] in [str(x) for x in range(101)]:

                    for _ in range(int(formula[x+1])):
                        y = list(i)
                        structure.extend(y)

        types = ["C", "H", "O", "N", "F", "S", "I", "Cl"]
        nums = {"S":0, "F":0, "C":0, "H": 0, "O" : 0, "N": 0, "I" : 0, "Cl" : 0}

        for i in formula:

            if i.upper() in types:
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

def compound():
    type_chem = input(("\n"*100)+"Enter the type of chemical compound (covalent or ionic): ").lower()
    init = input(("\n"*100)+"Enter the name of the chemical structure you wish to use: ").lower()
    info = getAll(init)
    chem_formula = info[0]
    real_structure = chem_formula[1]
    chem_formula = chem_formula[0]

    if type_chem == "covalent" or type_chem == "molecular":
        ans = ((createStructure(chem_formula, init).molecular()))
        predicted_structure, real_structure = str(ans).replace("\n", ""), real_structure.decode("utf8")
        print("\nThe predicted structure is: {}.\nThe real structure is {}".format(predicted_structure, real_structure))
        print("The formula is: " + str(chem_formula.decode("utf-8").replace("\n", "")))
        print("The png structure is opened.")
        time.sleep(5)

    elif type_chem == "ionic":
        ans = ((createStructure(chem_formula, init).ionic()))
        struct, charge = ans[0], ans[1]
        predicted_structure, real_structure = str(struct).replace("\n", ""), real_structure.decode("utf8")
        print("\nThe predicted structure is: {}.\nThe real structure is {}".format(predicted_structure, real_structure))
        print("The formula is: " + str(chem_formula.decode("utf-8").replace("\n", "")))
        print("The charge is: " + str(charge))
        print("The png structure is opened.")
        time.sleep(5)

def printed(stdscr, current_idx, typer, current_idy):
    stdscr.clear()
    h, w = stdscr.getmaxyx()
    counter = 0
    
    for idx, row in enumerate(typer):
        y = 0 +(idx % h)
        x = 0 + (math.floor(idx/h))*6
        if (x, y) == (current_idx, current_idy):
            stdscr.attron(curses.color_pair(1))
            stdscr.addstr(y, x, str(row))
            stdscr.attroff(curses.color_pair(1))
        else:
            stdscr.addstr(y, x, str(row))
    stdscr.refresh()

def periodic(stdscr):

    curses.curs_set(0)
    h, w = stdscr.getmaxyx()
    curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)

    current_row = 0
    current_pillar = 0
    printed(stdscr, current_row, menu, current_pillar)

    while True:

        for i in elements:
            comp = list(element_dict[i])
            sym = comp[0]
            current_symb = (stdscr.instr(current_pillar, (current_row+1), 2)).decode("utf8")
            current_symb = current_symb.replace(" ", "")
            if str(sym) == current_symb:
                name = i
                electconfig = comp[3]
                mass = comp[2]
                number = comp[1]
                stdscr.addstr(round(h/2)-3, (((math.floor(len(menu)//h))*6) + 12), "Name: " + str(name))
                stdscr.addstr(round(h/2)-4, (((math.floor(len(menu)//h))*6) + 12), "Electron Configuration: " + str(electconfig))
                stdscr.addstr(round(h/2)-5, (((math.floor(len(menu)//h))*6) + 12), "Mass: " + str(mass))
                stdscr.addstr(round(h/2)-6, (((math.floor(len(menu)//h))*6) + 12), "Atomic Number: " + str(number))

        stdscr.addstr(round(h/2), (((math.floor(len(menu)//h))*6) + 12), "Navigate with the arrow keys, and Press Q to quit at anytime.")
        stdscr.addstr(round(h/2)-1, (((math.floor(len(menu)//h))*6) + 12), "If cursor goes offscreen, hit the opposite arrow key to reverse.")
        
        c = stdscr.getch()
        if c == ord("q"):
            break

        elif c == curses.KEY_UP and current_pillar > 0:
            current_pillar -= 1
        elif c == curses.KEY_DOWN and current_pillar < h - 1:
            current_pillar += 1
        
        elif c == curses.KEY_LEFT and current_row > 0:
            current_row -= 6
        
        elif c == curses.KEY_RIGHT and current_row < (math.floor(len(menu)//h))*6:
            current_row += 6

        printed(stdscr, current_row, menu, current_pillar)
        stdscr.refresh()

    curses.endwin()


def stoich():
    pass




#main loop
while True:
    print("\n"*75)
    print(
            """
            
  _     _        ____ _                    
 | |__ (_) __ _ / ___| |__   ___ _ __ ___  
 | '_ \| |/ _` | |   | '_ \ / _ | '_ ` _ \ 
 | |_) | | (_| | |___| | | |  __| | | | | |
 |_.__/|_|\__, |\____|_| |_|\___|_| |_| |_|
          |___/                            
        """
)
    print("Compounds\nPeriodic Table\nStoichiometry (Coming soon!)")

    first = input("\nEnter what you want to use (type QUIT to exit): ")

    if first.lower() == "compounds" or first.lower() == "compound":
        compound()

    elif first.upper() == "QUIT":
        sys.exit()

    elif first.lower() == "periodic" or first.lower() == "periodic table":
        curses.wrapper(periodic)

