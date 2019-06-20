"""
BigC
bigChem - tools for chemistry
6.14.2019
Python3
"""


import requests
import json
import string
from io import BytesIO
from PIL import ImageTk, Image


def getAll(formula_name):
    url_formula = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/property/MolecularFormula/TXT"
    url_structure = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/property/CanonicalSMILES/TXT"
    url_png = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/PNG"
    response = requests.get(url=url_png)
    data_png = response.content
    img = (Image.open(BytesIO(data)))
    response = requests.get(url=url_formula)
    data_formula = response.content
    response = requests.get(url=url_structure)
    data_structure = response.content
    return list(data_formula, data_structure), img.show()


class createStructure:

    def __init__(self, formula, rname):
        self.formula = formula
        self.rname = rname
        formula = list(formula)

        for i in formula:
            if i == "\n":
                formula.remove(i)
            i = i.strip(" ")
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











init = input("Enter the name of the chemical structure you wish to use: ").lower()
chem_formula, real_structure = getAll(init)[0], getAll(init)[1]

ans = list((createStructure(chem_formula, init).molecular()))
predicted_structure, real_structure = str(ans[0]).replace("\n", ""), str(ans[1]).replace("\n", "")



print("\nThe predicted structure is: {}.\nThe real structure is {}.".format(predicted_structure, real_structure))
print("The formula is: " + str(chem_formula))
print("The png structure is opened.")
