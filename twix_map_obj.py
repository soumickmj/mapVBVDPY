import numpy as np
from copy import deepcopy

class twix_map_obj:
    """class to hold information about raw data from siemens MRI scanners
    (currently VB and VD software versions are supported and tested).
    """
    def __init__(self, fname, version, arg=None, dataType=None, rstraj=None):


        # self.sqzSize = None
        # self.sqzDims = None
        #
        # self.dataDims = None
        # self.NCol = None  # mdh information
        # self.NCha = None  # mdh information
        # self.NLin = None  # mdh information
        # self.NPar = None  # mdh information
        # self.NSli = None  # mdh information
        # self.NAve = None  # mdh information
        # self.NPhs = None  # mdh information
        # self.NEco = None  # mdh information
        # self.NRep = None  # mdh information
        # self.NSet = None  # mdh information
        # self.NSeg = None  # mdh information
        # self.NIda = None  # mdh information
        # self.NIdb = None  # mdh information
        # self.NIdc = None  # mdh information
        # self.NIdd = None  # mdh information
        # self.NIde = None  # mdh information
        # self.NAcq = None  # simple counter
        #
        # # mdh information
        # self.Lin = None
        # self.Par = None
        # self.Sli = None
        # self.Ave = None
        # self.Phs = None
        # self.Eco = None
        # self.Rep = None
        # self.Set = None
        # self.Seg = None
        # self.Ida = None
        # self.Idb = None
        # self.Idc = None
        # self.Idd = None
        # self.Ide = None
        #
        # self.centerCol = None
        # self.centerLin = None
        # self.centerPar = None
        # self.cutOff = None
        # self.coilSelect = None
        # self.ROoffcenter = None
        # self.timeSinceRF = None
        # self.IsReflected = None
        # self.IsRawDataCorrect = None  # SRY: storage for MDH flag raw data correct
        #
        # self.slicePos = None
        # self.freeParam = None
        # self.iceParam = None
        # self.scancounter = None
        # self.timestamp = None
        # self.pmutime = None
        # # memory position in file
        # self.memPos = None
        # self.isBrokenFile = None  # errors when parsing?
        #
        # self.arg = None  # arguments
        #
        # self.freadInfo = None
        # self.skipLin = None
        # self.skipPar = None

        self.dataType = 'image' if dataType is None else dataType.lower()
        self.filename = fname
        self.softwareVersion = version

        self.IsReflected = []
        self.IsRawDataCorrect = []
        self.NAcq = 0
        self.isBrokenFile = False

        self.dataDims = ['Col', 'Cha', 'Lin', 'Par', 'Sli', 'Ave', 'Phs', 'Eco', 'Rep', 'Set', 'Seg', 'Ida', 'Idb',
                         'Idc', 'Idd', 'Ide']
        self.freadInfo = {}
        self.arg = {}
        self.setDefaultFlags()
        if arg != None:
            self.arg |= arg

        self.flagAverageDim[self.dataDims.index('Ave')] = self.arg['doAverage']
        self.flagAverageDim[self.dataDims.index('Rep')] = self.arg['averageReps']
        self.flagAverageDim[self.dataDims.index('Set')] = self.arg['averageSets']
        self.flagAverageDim[self.dataDims.index('Seg')] = self.arg['ignoreSeg']

        if self.softwareVersion == 'vb':
            # every channel has its own full mdh
            self.freadInfo['szScanHeader'] = 0 # [bytes]
            self.freadInfo['szChannelHeader'] = 128 # [bytes]
            self.freadInfo['iceParamSz'] = 4

        elif self.softwareVersion == 'vd':
            if self.arg['doRawDataCorrect']:
                raise Exception('raw data correction for VD not supported/tested yet')

            self.freadInfo['szScanHeader'] = 192  # [bytes]
            self.freadInfo['szChannelHeader'] = 32  # [bytes]
            self.freadInfo['iceParamSz'] = 24 # vd version supports up to 24 ice params
        else:
            raise Exception('software version not supported')

        if rstraj is not None:
            self.rampSampTrj = rstraj
        else:
            self.rampSampTrj = []
            self.arg['rampSampRegrid'] = False

        # self.readerVersion
        # # flags
        # self.flagRemoveOS      # removes oversampling in read (col) during read operation
        # self.flagDoAverage         # averages over all avg during read operation
        # self.flagAverageReps             # averages over all repetitions
        # self.flagAverageSets            # averages over all sets
        # self.flagIgnoreSeg              # sum over all segments during read operation
        # self.flagSkipToFirstLine       # skips lines/partitions up to the first actually acquired line/partition
        #                                         # (e.g. only the center k-space is acquired in refscans, we don't want all the leading zeros in our data)
        #                                         # this is the default behaviour for everything but image scans (but can be changed manually)
        # self.flagRampSampRegrid         # perform on-the-fly ramp sampling regridding
        # self.flagDoRawDataCorrect        # SRY: apply raw data correction factors during read operation
        #
        # self.RawDataCorrectionFactors    # SRY: allow the user to set/get the factors
        #
        # self.flagAverageDim             # new: flags that determines whether certain dim. should be averaged/ignored
        # self.filename
        # self.softwareVersion
        # self.dataType
        # self.rampSampTrj
        #
        # #self.fullSize  # this is the full size of the data set according to the mdhs, i.e. flags
        # # like 'reduceOS' have no influence on it
        #
        # self.dataSize  # this is the current output size, depends on fullSize + some flags



    def __call__(self, *args, **kwargs):
        shape = kwargs.get('shape')
        bSqueeze = kwargs.get('squeeze', False)
        selRange, selRangeSz, outSize = self.calcRange(bSqueeze, shape=shape)

        # calculate page table (virtual to physical addresses)
        # this is now done every time, i.e. result is no longer saved in
        # a property - slower but safer (and easier to keep track of updates)
        ixToRaw, ixToTarget = self.calcIndices()
        tmp = np.reshape(list(range(np.product(self.fullSize[2:]).astype(np.int))), np.array(self.fullSize[2:], dtype=np.int), order='F')
        tmp = tmp.flatten(order='F')
        tmp = tmp[:selRangeSz[2:].prod()]
        tmp = tmp.reshape(selRangeSz[2:])
        # tmp = tmp[tuple(selRange[2:])] # Doesn't always work, so instead replaced with above 3 steps
        # tmp = tmp.squeeze()
        ixToRaw = ixToRaw[tmp]

        # delete all entries that point to zero (the "NULL"-pointer)
        notAcquired = (ixToRaw == 0)
        ixToRaw = ixToRaw[~notAcquired]
        ixToRaw = ixToRaw.astype(np.int)

        # calculate ixToTarg for possibly smaller, shifted + segmented
        # target matrix:
        cIx = np.ones((14, len(ixToRaw)))
        if ~self.flagAverageDim[2]:
            cIx[0, :] = self.Lin[ixToRaw] - self.skipLin
        if ~self.flagAverageDim[3]:
            cIx[1, :] = self.Par[ixToRaw] - self.skipPar
        if ~self.flagAverageDim[4]:
            cIx[2, :] = self.Sli[ixToRaw]
        if ~self.flagAverageDim[5]:
            cIx[3, :] = self.Ave[ixToRaw]
        if ~self.flagAverageDim[6]:
            cIx[4, :] = self.Phs[ixToRaw]
        if ~self.flagAverageDim[7]:
            cIx[5, :] = self.Eco[ixToRaw]
        if ~self.flagAverageDim[8]:
            cIx[6, :] = self.Rep[ixToRaw]
        if ~self.flagAverageDim[9]:
            cIx[7, :] = self.Set[ixToRaw]
        if ~self.flagAverageDim[10]:
            cIx[8, :] = self.Seg[ixToRaw]
        if ~self.flagAverageDim[11]:
            cIx[9, :] = self.Ida[ixToRaw]
        if ~self.flagAverageDim[12]:
            cIx[10, :] = self.Idb[ixToRaw]
        if ~self.flagAverageDim[13]:
            cIx[11, :] = self.Idc[ixToRaw]
        if ~self.flagAverageDim[14]:
            cIx[12, :] = self.Idd[ixToRaw]
        if ~self.flagAverageDim[15]:
            cIx[13, :] = self.Ide[ixToRaw]

        # make sure that indices fit inside selection range
        for k in range(2, len(selRange)):
            tmp = cIx[k - 2,:]
            for L in range(len(selRange[k])):
                cIx[k - 2, tmp==selRange[k][L]] = L

        sz = selRangeSz[2:]
        ixToTarg = self.sub2ind_double(sz, cIx[0, :], cIx[1, :], cIx[2, :],
                                       cIx[3, :], cIx[4, :], cIx[5, :], cIx[6, :],
                                       cIx[7, :], cIx[8, :], cIx[9, :], cIx[10, :],
                                       cIx[11, :], cIx[12, :], cIx[13, :])

        mem = self.memPos[ixToRaw]
        # sort mem for quicker access, sort cIxToTarg/Raw accordingly
        ix = np.argsort(mem)
        mem = np.sort(mem)
        ixToTarg = ixToTarg[ix]
        ixToRaw = ixToRaw[ix]

        return self.readData(mem,ixToTarg,ixToRaw,selRange,selRangeSz,outSize)

    def calcIndices(self):
        # calculate indices to target & source(raw)
        LinIx = self.Lin - self.skipLin
        ParIx = self.Par - self.skipPar
        sz = self.fullSize[2:]
        ixToTarget = self.sub2ind_double(sz, LinIx, ParIx,
                                         self.Sli, self.Ave, self.Phs,
                                         self.Eco, self.Rep, self.Set,
                                         self.Seg, self.Ida, self.Idb,
                                         self.Idc, self.Idd, self.Ide)

        # now calc. inverse index (page table: virtual to physical addresses)
        # indices of lines that are not measured are zero
        ixToRaw = np.zeros(int(np.product(self.fullSize[2:])))
        ixToRaw[ixToTarget] = list(range(len(ixToTarget)))
        return ixToRaw, ixToTarget

    def sub2ind_double(self, sz, *args):
        # SUB2IND_double Linear index from multiple subscripts
        ndx = args[-1] - 1
        for i in range(len(sz)-2, -1, -1):
            ix = args[i]
            ndx = (sz[i] * ndx) + ix - 1
        return ndx.astype(np.int)

    def calcRange(self, bSqueeze=False, shape=None):
        selRange = [[0] for _ in range((len(self.dataSize)))]
        outSize = np.ones(len(self.dataSize))

        if shape is None:
            for k in range(len(self.dataSize)):
                selRange[k] = list(range(self.dataSize[k]))
            outSize = self.sqzSize if bSqueeze else self.dataSize
        else:
            for k in range(len(shape)):
                cDim = self.dataDims.index(self.sqzDims[k]) if bSqueeze else k
                if shape[k] == -1:
                    if k < len(shape) - 1:
                        selRange[cDim] = list(range(int(self.dataSize[cDim])))
                    else: #  all later dimensions selected and 'vectorized'!
                        for L in range(cDim, len(self.dataSize)):
                            selRange[L] = list(range(int(self.dataSize[L])))

                        outSize[k] = np.prod(self.dataSize[cDim:])
                        break # jump out ouf for-loop
                elif isinstance(shape[k], int):
                    selRange[cDim] = [shape[k]-1]
                else:
                    raise ValueError(f'Unknown type {type(shape[k])} in shape')
                outSize[k] = len(selRange[cDim]) if isinstance(selRange[cDim], list) else 1
            for k in range(len(selRange)):
                if np.array(selRange[k]).max() > self.dataSize[k]:
                    raise ValueError('selection out of range')

        selRangeSz = np.ones(len(self.dataSize))
        for k in range(len(selRange)):
            selRangeSz[k] = len(selRange[k])

        # now select all indices for the dims that are averaged
        self.flagAverageDim = np.array(self.flagAverageDim, dtype=np.bool)
        for k in np.nonzero(self.flagAverageDim)[0]:
            selRange[k] = list(range(self.fullSize[k]))

        outSize = outSize.astype(np.int)
        selRangeSz = selRangeSz.astype(np.int)
        return selRange, selRangeSz, outSize

    def copy(self):
        """Copy function - to create object copies"""
        return deepcopy(self)

    def setDefaultFlags(self):
        """method to set flags to default values"""
        self.arg['removeOS'] = False
        self.arg['rampSampRegrid'] = False
        self.arg['doAverage'] = False
        self.arg['averageReps'] = False
        self.arg['averageSets'] = False
        self.arg['ignoreSeg'] = False
        self.arg['doRawDataCorrect'] = False
        self.flagAverageDim = [False] * 16

        self.arg['skipToFirstLine'] = self.dataType not in [
            'image',
            'phasecor',
            'phasestab',
        ]

        if 'rawDataCorrectionFactors' not in self.arg.keys():
            self.arg['rawDataCorrectionFactors'] = []

    def readMDH(self, mdh, filePos):
        """
        extract all values in all MDHs at once
        data types:
        Use double for everything non-logical, both ints and floats. Seems the
        most robust way to avoid unexpected cast-issues with very nasty side effects.
        Examples: eps(single(16777216)) == 2
                  uint32( 10 ) - uint32( 20 ) == 0
                  uint16(100) + 1e5 == 65535
                  size(1 : 10000 * uint16(1000)) ==  [1  65535]
        The 1st example always hits the timestamps.
        """

        if len(mdh) == 0 or not isinstance(mdh, dict):
            return

        self.NAcq = filePos.size
        sLC = np.double(mdh['sLC']) + 1  # +1: convert to matlab index style
        evalInfoMask1 = np.uint32(mdh['aulEvalInfoMask'][:, 0])

        # save mdh information for each line
        self.NCol = np.double(mdh['ushSamplesInScan'])
        self.NCha = np.double(mdh['ushUsedChannels'])
        self.Lin = sLC[:, 0]
        self.Ave = sLC[:, 1]
        self.Sli = sLC[:, 2]
        self.Par = sLC[:, 3]
        self.Eco = sLC[:, 4]
        self.Phs = sLC[:, 5]
        self.Rep = sLC[:, 6]
        self.Set = sLC[:, 7]
        self.Seg = sLC[:, 8]
        self.Ida = sLC[:, 9]
        self.Idb = sLC[:, 10]
        self.Idc = sLC[:, 11]
        self.Idd = sLC[:, 12]
        self.Ide = sLC[:, 13]

        self.centerCol = np.double(mdh['ushKSpaceCentreColumn']) + 1
        self.centerLin = np.double(mdh['ushKSpaceCentreLineNo']) + 1
        self.centerPar = np.double(mdh['ushKSpaceCentrePartitionNo']) + 1
        self.cutOff = np.double(mdh['sCutOff'])
        self.coilSelect = np.double(mdh['ushCoilSelect'])
        self.ROoffcenter = np.double(mdh['fReadOutOffcentre'])
        self.timeSinceRF = np.double(mdh['ulTimeSinceLastRF'])

        self.IsReflected = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 24), 1).astype(np.bool)

        self.scancounter = np.double(mdh['ulScanCounter'])
        self.timestamp = np.double(mdh['ulTimeStamp'])
        self.pmutime = np.double(mdh['ulPMUTimeStamp'])

        self.IsRawDataCorrect = np.minimum(np.bitwise_and(evalInfoMask1, 2 ** 10), 1) # SRY

        self.slicePos = np.double(mdh['SlicePos'])
        self.iceParam = np.double(mdh['aushIceProgramPara'])
        self.freeParam = np.double(mdh['aushFreePara'])

        self.memPos = filePos

    def tryAndFixLastMdh(self):
        eofWarning = [f'{__name__} :UnxpctdEOF']
        raise Warning(f'off {eofWarning}')

    def clean(self):
        if self.NAcq == 0:
            return

        # Cut mdh data to actual size. Maybe we rejected acquisitions at the end
        # due to read errors.

        fields = {'NCol', 'NCha',
                  'Lin', 'Par', 'Sli', 'Ave', 'Phs', 'Eco', 'Rep',
                  'Set', 'Seg', 'Ida', 'Idb', 'Idc', 'Idd', 'Ide',
                  'centerCol', 'centerLin', 'centerPar', 'cutOff',
                  'coilSelect', 'ROoffcenter', 'timeSinceRF', 'IsReflected',
                  'scancounter', 'timestamp', 'pmutime', 'IsRawDataCorrect',
                  'slicePos', 'iceParam', 'freeParam', 'memPos'}

        nack = self.NAcq
        idx = list(range(1, nack))

        for f in fields:
            if isinstance(self.__getattribute__(f), np.float):
                if nack < 1:
                    self.__setattr__(f, self.__getattribute__(f)[:idx])
            elif self.__getattribute__(f).shape[0] > nack: # rarely
                self.__setattr__(f, self.__getattribute__(f)[:idx]) # 1st dim: samples,  2nd dim acquisitions

        self.NLin = max(self.Lin)
        self.NPar = max(self.Par)
        self.NSli = max(self.Sli)
        self.NAve = max(self.Ave)
        self.NPhs = max(self.Phs)
        self.NEco = max(self.Eco)
        self.NRep = max(self.Rep)
        self.NSet = max(self.Set)
        self.NSeg = max(self.Seg)
        self.NIda = max(self.Ida)
        self.NIdb = max(self.Idb)
        self.NIdc = max(self.Idc)
        self.NIdd = max(self.Idd)
        self.NIde = max(self.Ide)

        # ok, let us assume for now that all NCol and NCha entries are
        # the same for all mdhs:
        if not isinstance(self.NCol, np.float):
            self.NCol = self.NCol[0]
        if not isinstance(self.NCha, np.float):
            self.NCha = self.NCha[0]

        if self.dataType == 'refscan' and self.NLin > 65500:
            self.Lin = ((self.Lin + (65536 - min(self.Lin[self.Lin > 65500], 0))) % 65536) + 1
            self.NLin = max(self.Lin)

        # to reduce the matrix sizes of non-image scans, the size
        # of the refscan_obj()-matrix is reduced to the area of the
        # actually scanned acs lines (the outer part of k-space
        # that is not scanned is not filled with zeros)
        # this behaviour is controlled by flagSkipToFirstLine which is
        # set to true by default for everything but image scans

        if not self.flagSkipToFirstLine:
            # the output matrix should include all leading zeros
            self.skipLin = 0
            self.skipPar = 0
        else:
            # otherwise, cut the matrix size to the start of the
            # first actually scanned line/partition (e.g. the acs/
            # phasecor data is only acquired in the k-space center)
            self.skipLin = min(self.Lin) - 1
            self.skipPar = min(self.Par) - 1

        NLinAlloc = max(1, self.NLin - self.skipLin)
        NParAlloc = max(1, self.NPar - self.skipPar)

        self.fullSize = np.array([self.NCol, self.NCha, NLinAlloc, NParAlloc,
                         self.NSli, self.NAve, self.NPhs, self.NEco,
                         self.NRep, self.NSet, self.NSeg, self.NIda,
                         self.NIdb, self.NIdc, self.NIdd, self.NIde], dtype=np.int)

        nByte = self.NCha * (self.freadInfo['szChannelHeader'] + (8 * self.NCol))

        # size for fread
        self.freadInfo['sz'] = [2, int(nByte / 8)]
        # reshape size
        self.freadInfo['shape'] = [int(self.NCol + (self.freadInfo['szChannelHeader'] / 8)), int(self.NCha)]
        # we need to cut MDHs from fread data
        self.freadInfo['cut'] = ((self.freadInfo['szChannelHeader'] / 8) + np.arange(0, self.NCol)).astype(np.int)

    def subsref(self, S):
        pass

    def unsorted(self, ival):
        pass

    def readData(self, mem, cIxToTarg=None, cIxToRaw=None, selRange=None, selRangeSz=None, outSize=None):

        if outSize is None:
            selRange[0] = -1
            selRange[1] = -1
            outSize = self.dataSize[0:2] + [len(mem)]
            selRangeSz = outSize
            cIxToTarg = range(0, selRangeSz[2])
            cIxToRaw = cIxToTarg
        else:
            if selRange[0] == list(range(0, int(self.dataSize[0]))):
                selRange[0] = -1
            if selRange[1] == list(range(0, int(self.dataSize[1]))):
                selRange[1] = -1
        out = np.zeros(outSize, dtype=np.complex)
        out = out.reshape(selRangeSz[0], selRangeSz[1], -1)

        if mem.size == 0:
            out = out.reshape(outSize)
            return

        cIxToTarg = self.cast2MinimalUint(cIxToTarg)

        szScanHeader = self.freadInfo['szScanHeader']
        readSize = self.freadInfo['sz']
        readShape = self.freadInfo['shape']
        readCut = self.freadInfo['cut']
        keepOS = list(range(int(self.NCol/4))) + list(range(int(self.NCol*3/4), int(self.NCol)))
        bRemoveOS = self.arg['removeOS']
        bIsReflected = self.IsReflected[cIxToRaw]
        bRegrid = self.flagRampSampRegrid and len(self.rampSampTrj)
        slicedata = self.slicePos[cIxToRaw, :].T

        # SRY store information about raw data correction
        bDoRawDataCorrect = self.arg['doRawDataCorrect']
        bIsRawDataCorrect = self.IsRawDataCorrect[cIxToRaw]
        isBrokenRead = False
        if bDoRawDataCorrect:
            rawDataCorrect = self.arg['rawDataCorrectionFactors']

        # MiV??: Raw data are read line-by-line in portions of 2xNColxNCha float32 points (2 for complex).
        # Computing and sorting(!) on these small portions is quite expensive, esp. when
        # it employs non-sequential memory paths. Examples are non-linear k-space acquisition
        # or reflected lines.
        # This can be sped up if slightly larger blocks of raw data are collected, first.
        # Whenever a block is full, we do all those operations and save it in the final "out" array.
        # What's a good block size? Depends on data size and machine (probably L2/L3/L4 cache sizes).
        # So...? Start with a small block, measure the time-per-line and double block size until
        # a minimum is found. Seems sufficiently robust to end up in a close-to-optimal size for every
        # machine and data.

        blockSz = 2 # size of blocks; must be 2^n; will be increased
        doLockblockSz = False # whether blockSZ should be left untouched
        tprev = float('inf') # previous time-per-line
        blockCtr = 0
        blockInit = np.full((int(readShape[0]), int(readShape[1]), blockSz), -np.inf)
        blockInit = blockInit.astype(np.complex)
        block = blockInit

        if bRegrid:
            v1 = list(range(selRangeSz[1]))
            v2 = list(range(blockSz))
            rsTrj = [self.rampSampTrj, v1, v2]
            trgTrj = np.linspace(start=min(self.rampSampTrj), stop=max(self.rampSampTrj), num=self.NCol)
            trgTrj = [trgTrj, v1, v2]

        # counter for proper scaling of averages/segments
        count_ave = np.zeros((1, 1, out.shape[2]))
        kMax = len(mem) # max loop index

        fid = self.fileopen()

        for k in range(kMax):
            # skip scan header
            fid.seek(int(mem[k] + szScanHeader), 0)

            raw = np.fromfile(fid, dtype=np.float32, count=int(np.product(readSize)))
            raw = raw.reshape(readSize, order='F').T

            # MiV??: With incomplete files fread() returns less than readSize points. The subsequent reshape will therefore error out.
            # We could check if len(raw) == np.prod(readSize), but people recommend exception handling for performance
            # reasons. Do it.
            try:
                raw_tmp = np.empty(raw.shape[:-1], dtype=np.complex)
                raw_tmp.real = raw[:,0]
                raw_tmp.imag = raw[:, 1]
                raw = raw_tmp.reshape(readShape, order='F')
            except:
                offset_bytes = mem[k] + szScanHeader
                raise Warning('An unexpected read error occurred at this byte offset: {0}. Will ignore this line and stop reading'.format(offset_bytes))

                # Reject this data fragment. To do so, init with the values of blockInit
                raw = np.full(shape=readShape, fill_value=-np.inf)
                isBrokenRead = true # remember it and bail out later

            block[:, :, blockCtr] = raw # fast serial storage in a cache array
            blockCtr = blockCtr + 1

            # Do expensive computations and reorderings on the gathered block.
            # Unfortunately, a lot of code is necessary, but that is executed much less
            # frequent, so its worthwhile for speed.
            # (was a MATLAB comment) TODO: Do *everything* block-by-block

            if (blockCtr == blockSz) or (k == kMax) or (isBrokenRead and blockCtr > 1):
                # remove MDH data from block:
                block = block[readCut,:,:]

                if bRegrid:
                    # correct for readout shifts
                    # the nco frequency is always scaled to the max.
                    # gradient amp and does account for ramp-sampling
                    ro_shift = self.calcOffcenterShiftRO(slicedata[:, k])
                    # TODO: come back here

                ix = list(range(1 + k - blockCtr, k + 1))

                if blockCtr != blockSz:
                    block = block[:, :, 0:blockCtr]

                if bRemoveOS: # remove oversampling in read
                    block = np.fft.ifft(block, axis=0)
                    block = np.fft.fft(block[keepOS, :, :], axis=0)

                if bDoRawDataCorrect and bIsRawDataCorrect[k]:
                    # SRY apply raw data correction if necessary
                    block = np.multiply(block, rawDataCorrect)

                isRefl = bIsReflected[ix]
                block[:, :, isRefl] = block[list(range(block.shape[0]-1, -1, -1)), :, :][:, :, isRefl]

                if (selRange[0] != -1) or (selRange[1] != -1):
                    block = block[selRange[0], selRange[1], :]

                I = np.argsort(cIxToTarg[ix])
                sortIdx = np.sort(cIxToTarg[ix])
                block = block[:, :, I] # reorder according to sorted target indices

                # Mark duplicate indices with 1; we'll have to treat them special for proper averaging
                # Bonus: The very first storage can be made much faster, because it's in-place.
                isDupe = np.array([False] + (np.diff(sortIdx)==0).tolist())

                idx1 = sortIdx[~isDupe] # acquired once in this block
                idxN = sortIdx[isDupe] # acquired multiple times

                count_ave[:, :, idx1] = count_ave[:, :, idx1] + 1

                if idxN.size == 0:
                    # no duplicates
                    if np.all(count_ave[:, :, idx1] == 1): # first acquisition of this line
                        out[:, :, idx1] = block

                    else:
                        out[:, :, idx1] = out[:, :, idx1] + block

                else:
                    out[:, :, idx1] = out[:, :, idx1] + block[:, :, ~isDupe]

                    block = block[:, :, isDupe]
                    for n in range(len(idxN)):
                        out[:, :, idxN[n]] = out[:, :, idxN[n]] + block[:, :, n]
                        count_ave[:, :, idxN[n]] = count_ave[:, :, idxN[n]] + 1

                # At the first few iterations, evaluate the spent time-per-line and decide
                # what to do with the block size.
                if ~doLockblockSz: # TODO: if speed problem -> refer to this portion of MATLAB code
                    # regression; reset size and lock it
                    blockSz = max(int(blockSz/2), 1)
                    blockInit = blockInit[:,:, 0:blockSz]
                    doLockblockSz = True

                    if bRegrid:
                        rsTrj[2] = list(range(blockSz))
                        trgTrj[2] = trgTrj[2]

                blockCtr = 0
                block = blockInit # reset to garbage

            if isBrokenRead:
                self.isBrokenFile = True
                break

        fid.close()

        # proper scaling (we don't want to sum our data but average it)
        # For large "out" bsxfun(@rdivide,out,count_ave) is incredibly faster than
        # bsxfun(@times,out,count_ave)!
        # @rdivide is also running in parallel, while @times is not. :-/

        if np.any(count_ave.reshape(-1, 1) > 1):
            count_ave = max(1, count_ave)
            out = np.divide(out, count_ave)

        out = out.reshape(outSize)
        out = out.squeeze()
        return out

    def calcOffcenterShiftRO(self, slicedata):
        # calculate ro offcenter shift from mdh's slicedata

        # slice position
        pos = slicedata[0:3]

        # quaternion
        a = slicedata[4]
        b = slicedata[5]
        c = slicedata[6]
        d = slicedata[3]

        read_dir = np.zeros((3,1))
        read_dir[0] = 2 * (a * b - c * d)
        read_dir[1] = 1 - 2 * (a^2 + c^2)
        read_dir[2] = 2 * (b * c + a * d)

        ro_shift = np.dot(pos, read_dir)

        return ro_shift

    def fileopen(self):
        fid = open(self.filename, mode='rb')
        return fid

    def cast2MinimalUint(self, N):
        Nmax = N.max()
        Nmin = N.min()

        if (Nmin < 0) or (Nmax > np.iinfo(np.uint64).max):
            return N

        if Nmax > np.iinfo(np.uint32).max:
            idxClass = np.uint64
        elif Nmax > np.iinfo(np.uint16).max:
            idxClass = np.uint32
        else:
            idxClass = np.uint16

        N = N.astype(idxClass)
        return N

    def resetFlags(self):
        pass

    @property
    def readerVersion(self):
        # returns utc-unixtime of last commit (from file precommit-unixtime)
        fid = open('precommit_unixtime', mode='rb')
        versiontime = int(fid.read())
        fid.close()
        return versiontime

    @property
    def flagRemoveOS(self):
        return self.arg['removeOS']

    @flagRemoveOS.setter
    def flagRemoveOS(self, val):
        # set method for removeOS
        self.arg['removeOS'] = bool(val)

    @property
    def flagDoAverage(self):
        ix = self.dataDims.index('Ave')
        return self.flagAverageDim[ix]

    @flagDoAverage.setter
    def flagDoAverage(self, val):
        ix = self.dataDims.index('Ave')
        self.flagAverageDim[ix] = val

    @property
    def flagAverageReps(self):
        ix = self.dataDims.index('Rep')
        return self.flagAverageDim[ix]

    @flagAverageReps.setter
    def flagAverageReps(self, val):
        ix = self.dataDims.index('Rep')
        self.flagAverageDim[ix] = val

    @property
    def flagAverageSets(self):
        ix = self.dataDims.index('Set')
        return self.flagAverageDim[ix]

    @flagAverageSets.setter
    def flagAverageSets(self, val):
        ix = self.dataDims.index('Set')
        self.flagAverageDim[ix] = val

    @property
    def flagIgnoreSeg(self):
        ix = self.dataDims.index('Seg')
        return self.flagAverageDim[ix]

    @flagIgnoreSeg.setter
    def flagIgnoreSeg(self, val):
        ix = self.dataDims.index('Seg')
        self.flagAverageDim[ix] = val

    @property
    def flagSkipToFirstLine(self):
        return self.arg['skipToFirstLine']

    @flagSkipToFirstLine.setter
    def flagSkipToFirstLine(self, val):
        val = bool(val)
        if val != self.arg['skipToFirstLine']:
            self.arg['skipToFirstLine'] = val

            if self.arg['skipToFirstLine']:
                self.skipLin = min(self.Lin) - 1
                self.skipPar = min(self.Par) - 1
            else:
                self.skipLin = 0
                self.skipPar = 0

            NLinAlloc = max(1, self.NLin - self.skipLin)
            NParAlloc = max(1, self.NPar - self.skipPar)
            self.fullSize[2:4] = [NLinAlloc, NParAlloc]

    @property
    def flagRampSampRegrid(self):
        return self.arg['rampSampRegrid']

    @flagRampSampRegrid.setter
    def flagRampSampRegrid(self, val):
        val = bool(val)
        if (val is True) and (len(self.rampSampTrj) == 0):
            raise Exception('No trajectory for regridding available')
        self.arg['rampSampRegrid'] = val

    @property
    def flagDoRawDataCorrect(self):
        return self.arg['doRawDataCorrect']

    @flagDoRawDataCorrect.setter
    def flagDoRawDataCorrect(self, val):
        val = bool(val)
        if (val is True) and (self.softwareVersion == 'vd'):
            raise Exception('raw data correction for VD not supported/tested yet')
        self.arg['doRawDataCorrect'] = val

    @property
    def RawDataCorrectionFactors(self):
        return self.arg['rawDataCorrectionFactors']

    @RawDataCorrectionFactors.setter
    def RawDataCorrectionFactors(self, val):
        # this may not work if trying to set the factors before NCha has
        # a meaningful value (ie before calling clean)
        if len(val) != self.NCha:
            raise Exception('RawDataCorrectionFactors must be a list of length NCha')

        self.arg['rawDataCorrectionFactors'] = val

    @property
    def dataSize(self):
        out = np.array(self.fullSize)

        if self.arg['removeOS']:
            ix = self.dataDims.index('Col')
            out[ix] = int(self.NCol/2)

        if self.flagAverageDim[0] or self.flagAverageDim[1]:
            raise Warning('averaging in col and cha dim not supported, resetting flag')
            self.flagAverageDim[0:2] = [False, False]

        out[self.flagAverageDim] = 1
        return out

    @property
    def sqzSize(self):
        return self.dataSize[self.dataSize > 1]

    @property
    def sqzDims(self):
        return np.array(self.dataDims)[self.dataSize > 1]

    @property
    def flagRampSampRegrid(self):
        return self.arg['rampSampRegrid']
