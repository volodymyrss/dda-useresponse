import glob
import os

import astropy.io.fits as fits
import numpy as np

import ast
import ddosa
import dataanalysis.core as da

import findic

class FindICARF(findic.FindICIndexEntry):
    ds="ISGR-ARF.-RSP"


class FindICRMF(findic.FindICIndexEntry):
    ds="ISGR-RMF.-RSP"

class FindResponse(ddosa.DataAnalysis):
#    input_findicarf = FindICARF
    input_findicrmf = FindICRMF

    def main(self):
        self.rmf_path=self.input_findicrmf.member_location

class FindICEBDS(findic.FindICIndexEntry):
    ds="ISGR-EBDS-MOD"

class SpectraBins(ddosa.SpectraBins):
    input_ic_ebds = FindICEBDS

    rmfbins=True

    version="v1"

    ebins=None

    def get_version(self):
        v=super(SpectraBins,self).get_version(self)
        if self.ebins is None:
            v += ".ic256"
        else:
            v += ".rebin%i"%len(self.ebins)

    def main(self):
        self.binrmf=self.input_ic_ebds.member_location
        e=fits.open(self.binrmf)[1].data
        self.bins=zip(e['E_MIN'],e['E_MAX'])
        self.binrmfext=self.binrmf+'[1]'

    def get_binrmfext(self):
        return self.binrmfext
