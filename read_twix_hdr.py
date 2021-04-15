import numpy as np
import regex as re

def read_twix_hdr(fid):
    # function to read raw data header information from siemens MRI scanners
    # (currently VB and VD software versions are supported and tested)
    nbuffers = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    rstraj = []

    prot = {}
    for b in range(nbuffers):
        # now read string up to null termination
        bufname = fid.read(10).decode(errors='ignore')
        bufname = re.findall('^\w*', bufname)
        bufname = bufname[0]
        fid.seek(fid.tell() + len(bufname) - 9, 0)
        buflen = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        buffer = fid.read(buflen).decode(errors='ignore')
        buffer = ''.join(re.split('\n\s*\n', buffer))
        prot[bufname] = parse_buffer(buffer)

    return prot, rstraj


def parse_buffer(buffer):
    ascconv = buffer.split('### ASCCONV BEGIN')
    if len(ascconv) > 1:
        ascconv = '\n'.join(ascconv[1].split('\n')[1:])
        ascconv = '\n'.join(ascconv.split('### ASCCONV END ###')[0].split('\n')[:-1])
        prot = parse_ascconv(ascconv)
        # problem in 5th iter call, COULD OR COULD NOT BE A PROBLEM
        xprot = buffer.split('### ASCCONV BEGIN')[0] + buffer.split('### ASCCONV END ###')[1][1:-1]

    else:
        ascconv = []
        prot = {}
        xprot = buffer

    if len(xprot) != 0:
        xprot = parse_xprot(xprot)
        prot.update(xprot)
    # xprot[list(prot.keys())[0]] = prot[list(prot.keys())[0]]
    # prot = xprot

    return prot
    # if not len(ascconv) == 0:
    #     ascconv


def parse_xprot(buffer):
    xprot = {}
    tokens = re.findall('<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)', buffer)
    tokens = tokens + re.findall('<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)', buffer)

    for m in range(len(tokens)):
        name = tokens[m][0]

        if not name[0].isalpha():
            name = 'x' + name

        value = tokens[m][-1]
        value = value.strip().replace('"', '')
        if value.isdecimal():
            value = int(value)
        elif value.isnumeric():
            value = float(value)

        xprot[name] = value

    return xprot


def parse_ascconv(buffer):
    mrprot = {}

    vararray = re.findall('(?P<name>\S*)\s*=\s*(?P<value>\S*)', buffer)
    tmp = {}
    for name, value in vararray:
        if value.isdecimal():
            value = int(value)
        elif value.isnumeric():
            value = float(value)

        if not '[' in name:
            continue

        p_name = name.split('.')[0].split('[')[0]
        if not p_name in tmp.keys():
            tmp[p_name] = {}
        ix = int(name.split('[')[1].split(']')[0])
        if ix not in tmp[p_name].keys():
            tmp[p_name][ix] = {}
        c_name = name.split('.')[-1]
        tmp[p_name][ix][c_name] = value

    mrprot[list(tmp.keys())[0]] = tmp[list(tmp.keys())[0]]

    return mrprot

        # v = re.findall('(?P<name>\w*)\[(?P<ix>[0-9]*)\]|(?P<name>\w*)', name)
        #
        # cnt = -1
        # tmp = {}
        # breaked = False
        # for var_name, ix in v:
        #     if not ((var_name == '') and (ix == '')):
        #         if not ix == '':
        #             if not var_name in tmp.keys():
        #                 tmp[var_name] = {}
        #             if not float(ix) in tmp[var_name].keys():
        #                 tmp[var_name][float(ix)] = {}
        #
        #
        #
        #
        #         if not var_name[0].isalpha():
        #             breaked = True
        #             break
        #
        #         cnt = cnt + 1
        #         tmp
        #         tmp.append({var_name})
        #
        #         if not ix == '':
        #             tmp[cnt][var_name] = []
        # if (not breaked) and (not len(tmp) == 0):
        #     S = tmp
