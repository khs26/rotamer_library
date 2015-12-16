import re

kcal_mol_in_joules = 4184
gas_constant = 8.3145

class RotamerMove():
    def __init__(self):
        """
        Initialise from either a file or a pickled RotamerMove instance.
        """
        self.residues = []
        self.deps = set()
        self.temp = None

    def read_settings(self, settings_filename):
        """
        Read settings for the rotamer move, which determines the residues that move and which distribution the move is
        selected from.

        :param settings_filename: filename of file containing settings
        """
        with open(settings_filename, "r") as settings_file:
            for line in settings_file:
                # List comp splits the string and gets rid of empty strings
                tokenized = [tok for tok in re.split("[ ,]", line) if tok]
                # Go through the possible options
                # Residue list
                if re.match("res", tokenized[0]):
                    self.residues = self.parse_res_list(" ".join(tokenized[1:]))
                # Dependencies
                elif any([re.match("depend", tok) for tok in tokenized]):
                    if any(["bb" in tok or "backbone" in tok for tok in tokenized]):
                        self.deps.add("backbone")
                    if any(["neighbour" in tok or "neighbor" in tok for tok in tokenized]):
                        self.deps.add("neighbour_idents")
                    if any(["sidechain" in tok or "sc" in tok for tok in tokenized]):
                        self.deps.add("neighbour_identities")
                        self.deps.add("neighbour_sidechains")
                # Rotamer distribution temperature
                elif re.match("temp", tokenized[0]):
                    kelvin = any(["K" in tok for tok in tokenized[1:]])
                    if kelvin:
                        tokenized[1:] = [tok.strip("K ") for tok in tokenized[1:] if tok.strip("K ")]
                        self.temp = float(tokenized[1]) * gas_constant / kcal_mol_in_joules
                    else:
                        self.temp = float(tokenized[1])

    def print_settings(self):
        """
        Prints the current settings for the rotamer moves.
        """
        for k, v in self.__dict__.items():
            print k, v

    def save_move(self):
        """
        Pickles the object, to save the move parameters between subsequent GMIN steps.

        :return:
        """

    def parse_res_list(self, res_list):
        """
        Parses the given residue list and returns a list of matching residue indices.

        :param res_list: list of residues according to the defined format
        :return: list of indices of matching residues
        """
        res_idxs = set()
        # Get rid of whitespace around hyphens, just in case
        res_list = re.sub(" - ", "-", res_list)
        # Split into tokens
        res_list = res_list.split()
        subtract = False
        for i, res in enumerate(res_list):
            if res == "not":
                subtract = True
                continue
            elif "-" in res:
                res_sublist = range(int(res.split("-")[0]), int(res.split("-")[1]) + 1)
            else:
                res_sublist = [int(res)]
            if subtract:
                for r in res_sublist:
                    res_idxs.remove(r)
                subtract = False
            else:
                for r in res_sublist:
                    res_idxs.add(r)
        return list(res_idxs)

    def __call__(self, coords):
        """
        Performs the rotamer move on a given set of coordinates.

        There are a few possible dependency levels:

        1) Identity of the amino acid being moved.

           In this case, we build a partition function by averaging over all of the neighbours and all neighbouring
           configurations.

        2) Identity of the amino acid and configuration of its backbone (phi/psi angles).

           In this case, we use variable kernel density estimates to construct conditional distributions of phi/psi
           angles as a function of rotamer states, and invert these using Bayes' rule to obtain P(r|phi, psi) (from
           Shapovalov MV and Dunbrack RL, Structure, 19, pp844-858, 2011)

        3) Identity of the amino acid and identity of the neighbouring amino acids.

           This is the same as case 1, but constructing separate partition functions for each pair of amino acid
           neighbours (i.e. 22 * 22 = 484 in total).

        4) Identity of the amino acid, its neighbours and rotameric states of the neighbours.

           This is the same as case 3, but constructing a conditional probability for the central rotamers in terms of
           the neighbouring rotameric states. *** TODO: Should this be an inversion of P(r_i-1, r_i+1 | r_i)? ***

        5) Combinations of 3 and 4 with backbone configurations.

        Each of these requires a slightly different method for calculating the probability distributions and different
        amounts of structural information prior to making a move.

        :param coords: Initial coordinates
        :return: New coordinates (after the move)
        """



if __name__ == "__main__":
    move = RotamerMove()
    move.read_settings("../tests/data/example_settings")
    move.print_settings()
