import os


class G09_Summary:
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_lines = open(self.file_path, "r").readlines()

        #  Some variables to save from output
        self.job_name = None
        self.atoms = None
        self.xyz = None
        self.theory = None
        self.output_summary = None
        self.vibrations_string = None
        self.HF = None
        self.MP2 = None
        self.Version = None
        self.ZPE = None
        self.Dipole = None
        self.DipoleDerivative = None
        self.NImag = None
        self.charge = None
        self.multiplicity = None
        self.num_atoms = None

    def record(self, line):
        self.output_summary.append(line)

    def saveVibrations(self, path=None):
        if self.vibrations_string is None:
            self.setVirbrations()

        if path:
            f_out = open(os.path.join(path, "{}_cart_disp.dat".format(self.job_name), "w"))
        else:
            f_out = open("{}_cart_disp.dat".format(self.job_name), "w")

        f_out.write(self.vibrations_string)
        f_out.close()

    def saveXYZ(self, path=None):
        if self.xyz is None:
            self.setInfo()

        if path:
            f_out = open(os.path.join(path, "{}.xyz".format(self.job_name), "w"))
        else:
            f_out = open("{}.xyz".format(self.job_name), "w")

        header = "{}\n\n".format(self.num_atoms)
        for coord in self.xyz:
            header += "{}    {}    {}    {}\n".format(*coord.split(","))

        f_out.write(header)
        f_out.close()

    def setVirbrations(self):
        start = False
        stop = False

        for i, line in enumerate(self.file_lines):
            if line.startswith("                      1                      2                      3"):
                start = i
            if start and line.__contains__("-------------------"):
                stop = i - 1
            if start and stop:
                break

        output_string = ""
        for x in self.file_lines[start:stop]:
            output_string += x
        self.vibrations_string = output_string

    def setInfo(self):
        self.output_summary = []
        reversed = list(self.file_lines)
        reversed.reverse()
        find_end = r"\\\@"
        find_start = " 1\\1\\"
        record = False

        for line in reversed:
            if line.startswith(find_start):
                record = False
            if line.__contains__(find_end):
                self.record(line)
                record = True
            if record:
                self.record(line)
        self.output_summary.reverse()

        output_string = ""
        for x in self.output_summary:
            output_string += x.rstrip()

        self.output_summary = output_string.rstrip()  # strips new line character

        output_string = ""
        for x in self.output_summary:
            output_string += x

        self.output_summary = output_string
        self.output_summary = self.output_summary.split("\\")
        self.output_summary = [i.replace(" ", "") for i in self.output_summary if i]

        self.job_name = self.output_summary[2]
        c_m = self.output_summary[3]
        self.charge = int(c_m.split(",")[0])
        self.multiplicity = int(c_m.split(",")[1])

        """
        Setting the XYZ coordinates
        """
        self.xyz = []
        for line in self.output_summary[4:]:
            if line.__contains__("Version"):
                break
            self.xyz.append(line)

        self.atoms = [x.split(",")[0].upper() for x in self.xyz]
        print(str(self.atoms)[1:-1])
        self.num_atoms = len(self.atoms)
        # for i, x in enumerate(self.output_summary):
        #     print(i, x)

    def generateNMS_F90(self):
        s = open("nms.f90.template", "r").read()
        atoms_formatted = str(self.atoms)[1:-1]
        f_out = open("{}_nms.f90".format(self.job_name), "w")
        f_out.write(s.format(self.num_atoms, atoms_formatted, "{}.xyz".format(self.job_name),
                             "{}_cart_disp.dat".format(self.job_name), self.job_name))
        f_out.close()


test = G09_Summary("butane_opt_freq.out")
test.saveXYZ()
test.saveVibrations()
test.generateNMS_F90()
