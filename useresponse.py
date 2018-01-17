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

    @property
    def ebds_mod_fn(self):
        return self.member_location

class CompressEBins(ddosa.DataAnalysis):
    input_ic_ebds=FindICEBDS

    factor=1

    def main(self):
        e = fits.open(self.input_ic_ebds.member_location)[1].data
        ic_bins = zip(e['E_MIN'], e['E_MAX'])

        print("raw bins",ic_bins)

        o_e1=e['E_MIN'][::int(self.factor)]
        o_e2=e['E_MAX'][self.factor-1::int(self.factor)]

        o_ebins=zip(o_e1,o_e2)
        o_e1,o_e2=zip(*o_ebins)

        print("compressed bins",o_ebins)

        f=fits.open(self.input_ic_ebds.member_location)
        f[1].data=f[1].data[:len(o_e1)]

        print(len(f[1].data['E_MIN']),len(o_e1))

        self.ebds_mod_fn = "compressed_ebins.fits"

        f[1].data['E_MIN'] = o_e1
        f[1].data['E_MAX'] = o_e2
        f.writeto(self.ebds_mod_fn,overwrite=True)




class SpectraBins(ddosa.SpectraBins):
    input_ic_ebds = CompressEBins

    rmfbins=True

    version="v1"

    def main(self):
        self.binrmf=self.input_ic_ebds.ebds_mod_fn
        e=fits.open(self.binrmf)[1].data
        self.bins=zip(e['E_MIN'],e['E_MAX'])
        self.binrmfext=self.binrmf+'[1]'

    def get_binrmfext(self):
        return self.binrmfext
