import glob
import os
import numpy as np
import datetime
from read_twix_hdr import read_twix_hdr
from twix_map_obj import twix_map_obj


class twix_obj_template:
    def __init__(self):
        self.hdr = None
        self.image = None
        self.noise = None
        self.phasecor = None
        self.phasestab = None
        self.phasestabRef0 = None
        self.phasestabRef1 = None
        self.refscan = None
        self.refscanPC = None
        self.refscanPS = None
        self.refscanPSRef0 = None
        self.refscanPSRef1 = None
        self.RTfeedback = None
        self.vop = None

def mapVBVD(filename, **kwargs):
    if filename == "":
        raise Exception('Please provide a valid filename')
        return
    else:
        if type(filename) == str:
            if not filename.endswith('.dat'):
                filename = filename + '.dat'
        else:
            # filename not a string, so assume that it is the MeasID
            measID = filename
            filelist = glob.glob('./*.dat')
            filesfound = 0
            for file in filelist:
                current_filename = file.split('\\')[-1]
                if 'MID' + str(measID).zfill(5) in current_filename:
                    if filesfound == 0:
                        filename = current_filename
                    filesfound = filesfound + 1
            if filesfound == 0:
                raise Exception('File with meas. id {0} not found.'.format(measID))
            elif filesfound > 1:
                print('Multiple files with meas. id {0} found. Choosing first occurence.'.format(measID))

    # add absolute path, when no path is given
    if os.path.dirname(filename) == '':
        filename = os.path.join(os.getcwd(), filename)

    # Parse varargin
    arg = {}
    arg['bReadImaScan'] = True
    arg['bReadNoiseScan'] = True
    arg['bReadPCScan'] = True
    arg['bReadRefScan'] = True
    arg['bReadRefPCScan'] = True
    arg['bReadRTfeedback'] = True
    arg['bReadPhaseStab'] = True
    arg['bReadHeader'] = True

    k = 1
    for key, value in kwargs.items():
        if not type(filename) == str:
            raise Exception('string expected')

        if key.lower() in ['readheader', 'readhdr', 'header', 'hdr']:
            arg['bReadHeader'] = value
        elif key.lower() in ['removeos', 'rmos']:
            arg['removeOS'] = value
        elif key.lower() in ['doaverage', 'doave', 'ave', 'average']:
            arg['doAverage'] = value
        elif key.lower() in ['averagereps', 'averagerepetitions']:
            arg['averageReps'] = value
        elif key.lower() in ['averagesets']:
            arg['averageSets'] = value
        elif key.lower() in ['ignseg', 'ignsegments', 'ignoreseg', 'ignoresegments']:
            arg['ignoreSeg'] = value
        elif key.lower() in ['rampsampregrid', 'regrid']:
            arg['rampSampRegrid'] = value
        elif key.lower() in ['rawdatacorrect', 'dorawdatacorrect']:
            arg['doRawDataCorrect'] = value
        else:
            raise Exception('Argument not recognized.')

    fid = open(filename, mode='rb')

    # get file size
    fid.seek(0, 2)
    fileSize = fid.tell()

    # start of actual measurement data (sans header)
    fid.seek(0, 0)

    firstInt = np.fromfile(fid, dtype=np.uint32, count=1)[0]
    secondInt = np.fromfile(fid, dtype=np.uint32, count=1)[0]

    # lazy software version check (VB or VD?)
    if (firstInt < 10000) and (secondInt <= 64):
        version = 'vd'
        print('Software version: VD (!?)')

        # number of different scans in file stored in 2nd in
        NScans = secondInt
        measID = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        fileID = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        measOffset = []
        measLength = []
        for k in range(NScans):
            measOffset.append(np.fromfile(fid, dtype=np.uint64, count=1)[0])
            measLength.append(np.fromfile(fid, dtype=np.uint64, count=1)[0])
            fid.seek(fid.tell() + 152 - 16, 0)
    else:
        # in VB versions, the first 4 bytes indicate the beginning of the raw data part of the file
        version = 'vb'
        print('Software version: VB (!?)')
        measOffset = [0]
        measLength = [fileSize]
        NScans = 1  # VB does not support multiple scans in one file

    # SRY read data correction factors
    # do this for all VB datasets, so that the factors are available later in the image_obj if the user chooses to set the correction flag
    if version is 'vb': # not implemented/tested for vd, yet
        datStart = measOffset[0] + firstInt
        fid.seek(0, 0)
        rawfactors = None
        while (fid.tell() < datStart) and (rawfactors == None):
            line = fid.readline()
            line = line.decode(errors='ignore')
            # line = np.frombuffer(line, dtype=np.uint32)

            # find the section of the protocol
            # Note: the factors are also available in <ParamArray."CoilSelects"> along with element name and FFT scale, but the parsing is significantly more difficult
            if '<ParamArray."axRawDataCorrectionFactor">' in line:
                while fid.tell() < datStart:
                    line = fid.readline()
                    line = line.decode(errors='ignore')
                    # find the line with correction factors, the factors are on the first line that begins with this pattern
                    if '{ {  { ' in line:
                        line = line.replace('}  { } ', '0.0')
                        line = line.replace('{', '')
                        line = line.replace('}', '')
                        rawfactors = np.loadtxt(line, dtype='float')
                        if (len(rawfactors) % 2) != 0:
                            raise Exception('Error reading rawfactors')
                        # note the transpose, this makes the vector multiplication during the read easier
                        arg['rawDataCorrectionFactors'] = complex(rawfactors[0:-1:2].T, rawfactors[1:-1:2].T)
                        break
        print('Read raw data correction factors')

    # data will be read in two steps (two while loops):
    # 1) reading all MDHs to find maximum line no., partition no.,... for ima, ref,... scan
    # 2) reading the data
    twix_obj = []

    for s in range(NScans):
        cPos = measOffset[s]
        fid.seek(cPos, 0)
        hdr_len = np.fromfile(fid, dtype=np.uint32, count=1)[0]

        # read header and calculate regridding (optional)
        rstraj = []
        if arg['bReadHeader']:
            twix_obj.append(twix_obj_template())
            twix_obj[s].hdr, rstraj = read_twix_hdr(fid)

        # declare data objects:
        twix_obj[s].image = twix_map_obj(arg = arg, dataType = 'image', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].noise= twix_map_obj(arg = arg, dataType = 'noise', fname = filename, version = version)
        twix_obj[s].phasecor = twix_map_obj(arg = arg, dataType = 'phasecor', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].phasestab = twix_map_obj(arg = arg, dataType = 'phasestab', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].phasestabRef0 = twix_map_obj(arg = arg, dataType = 'phasestab_ref0', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].phasestabRef1 = twix_map_obj(arg = arg, dataType = 'phasestab_ref1', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].refscan = twix_map_obj(arg = arg, dataType = 'refscan', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].refscanPC = twix_map_obj(arg = arg, dataType = 'refscan_phasecor', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].refscanPS = twix_map_obj(arg = arg, dataType = 'refscan_phasestab', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].refscanPSRef0 = twix_map_obj(arg = arg, dataType = 'refscan_phasestab_ref0', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].refscanPSRef1 = twix_map_obj(arg = arg, dataType = 'refscan_phasestab_ref1', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].RTfeedback = twix_map_obj(arg = arg, dataType = 'rtfeedback', fname = filename, version = version, rstraj = rstraj)
        twix_obj[s].vop = twix_map_obj(arg = arg, dataType = 'vop', fname = filename, version = version) # tx-array rf pulses

        if s == 0:
            # print('Reader version: {0}'.format(twix_obj[s]['image'].readerVersion))
            date_str = datetime.datetime.utcfromtimestamp(1496740353).ctime()
            print('(UTC: {0})'.format(date_str))

        # jump to first mdh
        cPos = cPos + hdr_len
        fid.seek(cPos, 0)

        # find all mdhs and save them in binary form, first:
        print('Scan {0}/{1}, read all mdhs:'.format(s+1, NScans))

        mdh_blob, filePos, isEOF = loop_mdh_read(fid, version, NScans, s, measOffset[s], measLength[s]) # uint8; size: [byteMDH  Nmeas]

        cPos = filePos[-1]
        filePos = filePos[:-1]

        # get mdhs and masks for each scan, no matter if noise, image, RTfeedback etc:
        mdh, mask = evalMDH(mdh_blob, version) # this is quasi-instant (< 1s) :-)

        # Assign mdhs to their respective scans and parse it in the correct twix objects.

        if arg['bReadImaScan']:
            tmpMdh = {}
            isCurrScan = mask['MDH_IMASCAN']
            isCurrScan = isCurrScan.flatten()
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].image.readMDH(tmpMdh, filePos[isCurrScan == 1])

        if arg['bReadNoiseScan']:
            tmpMdh = {}
            isCurrScan = mask['MDH_NOISEADJSCAN']
            isCurrScan = isCurrScan.flatten()
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].noise.readMDH(tmpMdh, filePos[isCurrScan == 1])

        if arg['bReadRefScan']:
            tmpMdh = {}
            isCurrScan = (mask['MDH_PATREFSCAN'] | mask['MDH_PATREFANDIMASCAN']) & \
                         (~(mask['MDH_PHASCOR'] | mask['MDH_PHASESTABSCAN'] |
                            mask['MDH_REFPHASESTABSCAN'] | mask['MDH_RTFEEDBACK'] | mask['MDH_HPFEEDBACK']))
            isCurrScan = isCurrScan.flatten()
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].refscan.readMDH(tmpMdh, filePos[isCurrScan == 1])

        if arg['bReadRTfeedback']:
            tmpMdh = {}
            isCurrScan = (mask['MDH_RTFEEDBACK'] | mask['MDH_HPFEEDBACK']) & (~(mask['MDH_VOP']))
            isCurrScan = isCurrScan.flatten()
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].RTfeedback.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = (mask['MDH_RTFEEDBACK'] & mask['MDH_VOP'])
            isCurrScan = isCurrScan.flatten()
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].vop.readMDH(tmpMdh, filePos[isCurrScan == 1])

        if arg['bReadPCScan']:
            # logic really correct?
            tmpMdh = {}
            isCurrScan = mask['MDH_PHASCOR'] & (~(mask['MDH_PATREFSCAN']) | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].phasecor.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = mask['MDH_PHASCOR'] & (mask['MDH_PATREFSCAN'] | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].refscanPC.readMDH(tmpMdh, filePos[isCurrScan == 1])

        if arg['bReadPhaseStab']:
            tmpMdh = {}
            isCurrScan = (mask['MDH_PHASESTABSCAN'] & (~mask['MDH_REFPHASESTABSCAN'])) & (~(mask['MDH_PATREFSCAN']) | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].phasestab.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = (mask['MDH_PHASESTABSCAN'] & (~mask['MDH_REFPHASESTABSCAN'])) & (mask['MDH_PATREFSCAN'] | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].refscanPS.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = (mask['MDH_REFPHASESTABSCAN'] & (~mask['MDH_PHASESTABSCAN'])) & (~(mask['MDH_PATREFSCAN']) | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].phasestabRef0.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = (mask['MDH_REFPHASESTABSCAN'] & (~mask['MDH_PHASESTABSCAN'])) & (mask['MDH_PATREFSCAN'] | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].refscanPSRef0.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = (mask['MDH_REFPHASESTABSCAN'] & mask['MDH_PHASESTABSCAN']) & (~(mask['MDH_PATREFSCAN']) | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].phasestabRef1.readMDH(tmpMdh, filePos[isCurrScan == 1])

            tmpMdh = {}
            isCurrScan = (mask['MDH_REFPHASESTABSCAN'] & mask['MDH_PHASESTABSCAN']) & (mask['MDH_PATREFSCAN'] | mask['MDH_PATREFANDIMASCAN'])
            for f in mdh.keys():
                tmpMdh[f] = mdh[f][isCurrScan == 1]
            twix_obj[s].refscanPSRef1.readMDH(tmpMdh, filePos[isCurrScan == 1])

        mdh = {}
        tmpMdh = {}
        filePos = []
        isCurrScan = []

        for scan in ['image', 'noise', 'phasecor', 'phasestab',
                     'phasestabRef0', 'phasestabRef1', 'refscan',
                     'refscanPC', 'refscanPS', 'refscanPSRef0',
                     'refscanPSRef1', 'RTfeedback', 'vop']:
            f = scan
            # remove unused fields
            if twix_obj[s].__getattribute__(f).NAcq == 0:
                twix_obj[s].__delattr__(f)
            else:
                if isEOF:
                    # recover from read error
                    twix_obj[s].__getattribute__(f).tryAndFixLastMdh()
                else:
                    twix_obj[s].__getattribute__(f).clean()

    # end of NScans loop

    if NScans == 1:
        twix_obj = twix_obj[0]

    return twix_obj
# end of mapVBVD()

def loop_mdh_read(fid, version, Nscans, scan, measOffset, measLength):
    # Goal of this function is to gather all mdhs in the dat file and store them
    # in binary form, first. This enables us to evaluate and parse the stuff in
    # a MATLAB-friendly (vectorized) way. We also yield a clear separation between
    # a lengthy loop and other expressions that are evaluated very few times.
    #
    # The main challenge is that we never know a priori, where the next mdh is
    # and how many there are. So we have to actually evaluate some mdh fields to
    # find the next one.
    #
    # All slow things of the parsing step are found in the while loop.
    # => It is the (only) place where micro-optimizations are worthwhile.
    #
    # The current state is that we are close to sequential disk I/O times.
    # More fancy improvements may be possible by using workers through parfeval()
    # or threads using a java class (probably faster + no toolbox):
    # http://undocumentedmatlab.com/blog/explicit-multi-threading-in-matlab-part1

    if version == 'vb':
        isVD = False
        byteMDH = 128
    elif version == 'vd':
        isVD = True
        byteMDH = 184
        szScanHeader = 192 # [bytes]
        szChannelHeader = 32 # [bytes]
    else:
        # arbitrary assumptions:
        isVD = False
        byteMDH = 128
        raise Warning('UnknownVer. Software version {0} is not supported.'.format(version))

    cPos = fid.tell()
    n_acq = 0
    allocSize = 4096
    ulDMALength = byteMDH
    isEOF = False
    last_progress = 0

    mdh_blob = np.zeros((byteMDH,0), dtype='uint8')
    szBlob = mdh_blob.shape[1]
    filePos = np.zeros((0, 1), dtype=type(cPos))

    fid.seek(cPos, 0)

    # ======================================
    # constants and conditional variables
    # ======================================
    bit_0 = 2 ** 0
    bit_5 = 2 ** 5
    mdhStart = - byteMDH

    u8_000 = np.zeros(3, dtype='uint8') # for comparison with data_u8(1:3)

    # 20 fill bytes in VD (21:40)
    evIdx = 21 + (20 * isVD) - 1 # 1st byte of evalInfoMask
    dmaIdx = [x + (20 * isVD) - 1 for x in range(29, 33)] # to correct DMA length using NCol and NCha

    if isVD:
        dmaOff = szScanHeader
        dmaSkip = szChannelHeader
    else:
        dmaOff = 0
        dmaSkip = byteMDH
    # ======================================

    while True:
        '''Read mdh as binary (uint8) and evaluate as little as possible to know
           where the next mdh is (ulDMALength / ushSamplesInScan & ushUsedChannels)
           whether it is only for sync (MDH_SYNCDATA)
           whether it is the last one (MDH_ACQEND)
           evalMDH() contains the correct and readable code for all mdh entries.'''

        try:
            # read everything and cut out the mdh
            data_u8 = np.fromfile(fid, dtype=np.uint8, count=int(ulDMALength))
            data_u8 = data_u8[mdhStart:]
        except:
            raise Warning('UnxpctdEOF. An unexpected read error occurred at this byte offset: {0} ({1} GiB). Will stop reading now'.format(cPos, cPos/1024**3))
            isEOF = True
            break

        bitMask = data_u8[evIdx] # the initial 8 bit from evalInfoMask are enough
        #print(bitMask)
        if (np.all(data_u8[0:3] == u8_000)) or (bitMask & bit_0):
            # ok, look closer if really all *4* bytes are 0:
            data_u8[4] = (data_u8[4] % 2) # ubit24: keep only 1 bit from the 4th byte
            ulDMALength = data_u8[0:4].view(np.uint32).astype('double')[0]

            if (ulDMALength == 0) or (bitMask & bit_0):
                cPos = cPos + ulDMALength
                # jump to next full 512 bytes
                if cPos % 512:
                    cPos = cPos + 512 - (cPos % 512)
                break
        if (bitMask & bit_5): # MDH_SYNCDATA
            data_u8[3] = (data_u8[3] % 2) # ubit24: keep only 1 bit from the 4th byte
            ulDMALength = data_u8[0:4].view(np.uint32).astype('double')[0]
            cPos = cPos + ulDMALength
            continue

        # pehses: the pack bit indicates that multiple ADC are packed into one
        # DMA, often in EPI scans (controlled by fRTSetReadoutPackaging in IDEA)
        # since this code assumes one adc (x NCha) per DMA, we have to correct
        # the "DMA length"
        #       if mdh.ulPackBit
        # it seems that the packbit is not always set correctly
        NCol_NCha = data_u8[dmaIdx].view(np.uint16).astype('double') # [ushSamplesInScan  ushUsedChannels]
        ulDMALength = dmaOff + (8 * NCol_NCha[0] + dmaSkip) * NCol_NCha[1]

        n_acq = n_acq + 1

        # grow arrays in batches
        if n_acq > szBlob:
            if mdh_blob.shape[1] == 0:
                mdh_blob = np.zeros((mdh_blob.shape[0], allocSize))
                filePos = np.zeros(allocSize)
            else:
                mdh_blob = np.concatenate((mdh_blob, np.zeros((mdh_blob.shape[0], allocSize))), axis=1)
                filePos = np.concatenate((filePos, np.zeros(allocSize)))
            szBlob = mdh_blob.shape[1]

        mdh_blob[:, n_acq-1] = data_u8
        filePos[n_acq-1] = cPos

        cPos = cPos + ulDMALength

    if isEOF:
        n_acq = n_acq - 1 # ignore the last attempt

    filePos[n_acq] = cPos # save pointer to the next scan

    # discard overallocation:
    mdh_blob = mdh_blob[:, 0:n_acq]
    filePos = filePos[0:n_acq+1].T # row vector

    print('{0} MB read'.format(measLength / (1024 ** 2)))

    mdh_blob = mdh_blob.astype('uint8')
    return mdh_blob, filePos, isEOF


def evalMDH(mdh_blob, version):
    # see pkg/MrServers/MrMeasSrv/SeqIF/MDH/mdh.h
    # and pkg/MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h
    mdh = {}
    mask = {}

    if not mdh_blob.dtype == np.uint8:
        raise Exception('mdh data must be a uint8 array!')

    if version[-1] == 'd':
        isVD = True
        # mdh_blob =np.delete(mdh_blob, range(21, 41), axis=0) # remove 20 unnecessary bytes
        mdh_blob = np.delete(mdh_blob, range(20, 40), axis=0)  # remove 20 unnecessary bytes
    else:
        isVD = False

    Nmeas = mdh_blob.shape[1]

    mdh['ulPackBit'] = (mdh_blob[3, :] >> 2) & 1
    mdh['ulPCI_rx'] = set_bit(set_bit(mdh_blob[3, :], 7, 0), 8, 0) # keep 6 relevant bits
    mdh_blob[3, :] = (mdh_blob[3, :] >> 1) & 1 # ubit24: keep only 1 bit from the 4th byte

    data_uint32 = mdh_blob[0:76, :].T.flatten().view(np.uint32)
    data_uint16 = mdh_blob[28:, :].T.flatten().view(np.uint16)
    data_single = mdh_blob[68:, :].T.flatten().view(np.single)

    data_uint32 = data_uint32.reshape(Nmeas, -1)
    data_uint16 = data_uint16.reshape(Nmeas, -1)
    data_single = data_single.reshape(Nmeas, -1)
                                                            # byte pos.
    # mdh['ulDMALength'] = data_uint32[:, 0]                # 1 :   4
    mdh['lMeasUID'] = data_uint32[:, 1]                     # 5 :   8
    mdh['ulScanCounter'] = data_uint32[:, 2]                # 9 :   12
    mdh['ulTimeStamp'] = data_uint32[:, 3]                  # 13 :  16
    mdh['ulPMUTimeStamp'] = data_uint32[:, 4]               # 17 :  20
    mdh['aulEvalInfoMask'] = data_uint32[:, 5:7]            # 21 :  28
    mdh['ushSamplesInScan'] = data_uint16[:, 0]             # 29 :  30
    mdh['ushUsedChannels'] = data_uint16[:, 1]              # 31 :  32
    mdh['sLC'] = data_uint16[:, 2:16]                       # 33 :  60
    mdh['sCutOff'] = data_uint16[:, 16:18]                  # 61 :  64
    mdh['ushKSpaceCentreColumn'] = data_uint16[:, 18]       # 66 :  66
    mdh['ushCoilSelect'] = data_uint16[:, 19]               # 67 :  68
    mdh['fReadOutOffcentre'] = data_single[:, 0]            # 69 :  72
    mdh['ulTimeSinceLastRF'] = data_uint32[:, 18]           # 73 :  76
    mdh['ushKSpaceCentreLineNo'] = data_uint16[:, 24]       # 77 :  78
    mdh['ushKSpaceCentrePartitionNo'] = data_uint16[:, 25]  # 79 :  80

    if isVD:
        mdh['SlicePos'] = data_single[:, 3:10]              # 81 : 108
        mdh['aushIceProgramPara'] = data_uint16[:, 40:64]   # 109 : 156
        mdh['aushFreePara'] = data_uint16[:, 64:68]         # 157 : 164
    else:
        mdh['aushIceProgramPara'] = data_uint16[:, 26:30]   # 81 :  88
        mdh['aushFreePara'] = data_uint16[:, 30:34]         # 89 :  96
        mdh['SlicePos'] = data_single[:, 7:14]              # 97 : 124

    # # inlining of evalInfoMask
    # evalInfoMask1 = mdh['aulEvalInfoMask'][:, 0]
    # mask['MDH_ACQEND'] = min(np.bitwise_and(evalInfoMask1, 2 ** 0), 1)
    # mask['MDH_RTFEEDBACK'] = min(np.bitwise_and(evalInfoMask1, 2 ** 1), 1)
    # mask['MDH_HPFEEDBACK'] = min(np.bitwise_and(evalInfoMask1, 2 ** 2), 1)
    # mask['MDH_SYNCDATA'] = min(np.bitwise_and(evalInfoMask1, 2 ** 5), 1)
    # mask['MDH_RAWDATACORRECTION'] = min(np.bitwise_and(evalInfoMask1, 2 ** 10), 1)
    # mask['MDH_REFPHASESTABSCAN'] = min(np.bitwise_and(evalInfoMask1, 2 ** 14), 1)
    # mask['MDH_PHASESTABSCAN'] = min(np.bitwise_and(evalInfoMask1, 2 ** 15), 1)
    # mask['MDH_SIGNREV'] = min(np.bitwise_and(evalInfoMask1, 2 ** 17), 1)
    # mask['MDH_PHASCOR'] = min(np.bitwise_and(evalInfoMask1, 2 ** 21), 1)
    # mask['MDH_PATREFSCAN'] = min(np.bitwise_and(evalInfoMask1, 2 ** 22), 1)
    # mask['MDH_PATREFANDIMASCAN'] = min(np.bitwise_and(evalInfoMask1, 2 ** 23), 1)
    # mask['MDH_REFLECT'] = min(np.bitwise_and(evalInfoMask1, 2 ** 24), 1)
    # mask['MDH_NOISEADJSCAN'] = min(np.bitwise_and(evalInfoMask1, 2 ** 25), 1)
    # mask['MDH_VOP'] = min(np.bitwise_and(mdh['aulEvalInfoMask'][2], 2 ** (53-32)), 1) # was 0 in VD
    # mask['MDH_IMASCAN'] = np.ones((Nmeas, 1), dtype=np.uint32)
    #
    # noImaScan = mask['MDH_ACQEND'] | mask['MDH_RTFEEDBACK'] | mask['MDH_HPFEEDBACK'] | \
    #             mask['MDH_PHASCOR'] | mask['MDH_NOISEADJSCAN'] | mask['MDH_PHASESTABSCAN'] | \
    #             mask['MDH_REFPHASESTABSCAN'] | mask['MDH_SYNCDATA'] | mask['MDH_PATREFSCAN'] & \
    #             (not mask['MDH_PATREFANDIMASCAN'])
    #
    # mask['MDH_IMASCAN'][noImaScan] = 0


    # inlining of evalInfoMask
    evalInfoMask1 = mdh['aulEvalInfoMask'][:, 0]
    mask['MDH_ACQEND'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 0), 1)
    mask['MDH_RTFEEDBACK'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 1), 1)
    mask['MDH_HPFEEDBACK'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 2), 1)
    mask['MDH_SYNCDATA'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 5), 1)
    mask['MDH_RAWDATACORRECTION'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 10), 1)
    mask['MDH_REFPHASESTABSCAN'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 14), 1)
    mask['MDH_PHASESTABSCAN'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 15), 1)
    mask['MDH_SIGNREV'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 17), 1)
    mask['MDH_PHASCOR'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 21), 1)
    mask['MDH_PATREFSCAN'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 22), 1)
    mask['MDH_PATREFANDIMASCAN'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 23), 1)
    mask['MDH_REFLECT'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 24), 1)
    mask['MDH_NOISEADJSCAN'] = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 25), 1)
    mask['MDH_VOP'] = np.minimum(np.bitwise_and(mdh['aulEvalInfoMask'][1][0], 2 ** (53-32)), 1) # was 0 in VD
    mask['MDH_IMASCAN'] = np.ones((Nmeas, 1), dtype=np.uint32)

    noImaScan = mask['MDH_ACQEND'] | mask['MDH_RTFEEDBACK'] | mask['MDH_HPFEEDBACK'] | \
                mask['MDH_PHASCOR'] | mask['MDH_NOISEADJSCAN'] | mask['MDH_PHASESTABSCAN'] | \
                mask['MDH_REFPHASESTABSCAN'] | mask['MDH_SYNCDATA'] | (mask['MDH_PATREFSCAN'] & (~mask['MDH_PATREFANDIMASCAN']))

    mask['MDH_IMASCAN'][noImaScan == 1] = 0


    return mdh, mask


# Utility function
def set_bit(v, index, x):
  """Set the index:th bit of v to 1 if x is truthy, else to 0, and return the new value."""
  mask = 1 << index   # Compute mask, an integer with just bit 'index' set.
  v & ~mask          # Clear the bit indicated by the mask (if x is False)
  if x:
    v | mask         # If x was True, set the bit indicated by the mask.
  return v



if __name__ == '__main__':
    twix_obj = mapVBVD('meas_MID00051_FID11835_benz_vibe_cap3_qfat_inhale.dat')
    # mapVBVD(51)
