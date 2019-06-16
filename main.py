"""
BigC
Generative Art w/ Chemical Structures
6.14.2019
Python3
"""


import requests
import json
import string
from io import BytesIO
from PIL import ImageTk, Image


def getFormula(formula_name):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + formula_name + "/property/MolecularFormula/TXT"
    response = requests.get(url=url)
    content = str(response.text)
    return content


def getPNG(bruh):
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/" + bruh + "/PNG"
    response = requests.get(url=url)
    data = response.content
    img = (Image.open(BytesIO(data)))
    return img.show()


class createStructure:

    def __init__(self, formula):
        self.formula = formula
        formula = list(formula)

        for i in formula:
            if i == "\n":
                formula.remove(i)
            i = i.strip(" ")
        self.formula = formula



    def molecular(self):

        formula = self.formula
        structure = []

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
                x = 0
                y = []

                if len(formula[:z+x+2]) == len(formula):
                    nums[i] = int(formula[-1])
                    break

                while formula[z+1+x] not in list(string.ascii_uppercase):


                    y.extend(formula[z+1+x:z+2+x])

                    x += 1



                if int("".join(y)) > 1:
                    nums[i] = int("".join(y))

                else:
                    nums[i] = 1


        count = 2

        for c, i in enumerate(structure):

            if int(nums["H"]) != 1 and int(nums["O"]) != 0 and count != 0:
                    structure[c] = list(str(i) + "HOH")
                    nums["H"] -= 2
                    nums["O"] -= 1
                    count -= 1

            elif int(nums["H"]) - 3 > -1 and count != 0:
                    structure[c] = list(i + ("H"*3))
                    nums["H"] -= 3

            elif nums["O"] != 0 and nums["H"] != 0:
                    structure[c] = list(i + ("OH"))
                    nums["O"] -= 1
                    nums["H"] -= 1

            elif nums["F"] - 3 > -1 and count != 0:
                structure[c] = list(i + ("F"*3))
                nums["F"] -= 3
                count -= 1

            elif nums["F"] - 2 > -1:
                structure[c] = list(i + ("F"*2))
                nums["F"] -= 2

            elif nums["H"] -2 > -1:
                structure[c] = list(i + ("H"*2))
                nums["H"] -= 2






        if len(structure) == 0:
            structure = formula

        return structure



    def checker(self):
        formula = self.formula
        if "C" in formula:
            createStructure(formula).molecular()









init = input("Enter the name of the chemical structure you wish to use: ").lower()
chem_formula = getFormula(init)
png = getPNG(init)

structure = (createStructure(chem_formula).molecular())
print("The structure is: " + str(structure))
print("The formula is: " + str(chem_formula))
print("The png structure is opened.")
