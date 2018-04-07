import glob
import os

import astropy.io.fits as fits
import numpy as np

import ast
import ddosa
import dataanalysis.core as da
import pilton

from dataanalysis.importing import load_by_name
import findic
#findic=load_by_name("git://findic/icversion")[0]

class FindICARF(findic.FindICIndexEntry):
    ds="ISGR-ARF.-RSP"


class FindICRMF(findic.FindICIndexEntry):
    ds="ISGR-RMF.-RSP"

class FindResponse(ddosa.DataAnalysis):
#    input_findicarf = FindICARF
    input_scw=ddosa.ScWData
    input_findicrmf = FindICRMF

    def main(self):
        self.rmf_path=self.input_findicrmf.get_member_location(self.input_scw)

class FindICEBDS(findic.FindICIndexEntry):
    ds="ISGR-EBDS-MOD"


class CompressEBins(ddosa.DataAnalysis):
    input_ic_ebds=FindICEBDS

    factor=1

    def get_version(self):
        v=super(CompressEBins,self).get_version()
        v+="_f%i"%self.factor
        return v

    def main(self):
        e = fits.open(self.input_ic_ebds.get_member_location())[1].data
        ic_bins = zip(e['E_MIN'], e['E_MAX'])

        print("raw bins",ic_bins)

        o_e1=e['E_MIN'][::int(self.factor)]
        o_e2=e['E_MAX'][self.factor-1::int(self.factor)]

        o_ebins=zip(o_e1,o_e2)
        o_e1,o_e2=zip(*o_ebins)

        print("compressed bins",o_ebins)

        f=fits.open(self.input_ic_ebds.get_member_location())
        f[1].data=f[1].data[:len(o_e1)]

        print(len(f[1].data['E_MIN']),len(o_e1))

        self.ebds_mod_fn = "compressed_ebins.fits"

        f[1].data['E_MIN'] = o_e1
        f[1].data['E_MAX'] = o_e2
        f.writeto(self.ebds_mod_fn,overwrite=True)


class RebinResponse(ddosa.DataAnalysis):
    input_rsp = FindResponse
    input_ic_ebds = FindICEBDS
    input_ebins=CompressEBins

    def main(self):
        new_e = fits.open(self.input_ebins.ebds_mod_fn)[1]
        new_bins = zip(new_e.data['E_MIN'], new_e.data['E_MAX'])
        print("new bins:",new_bins)

        orig_e = fits.open(self.input_ic_ebds.get_member_location())[1]
        orig_bins = zip(orig_e.data['E_MIN'], orig_e.data['E_MAX'])
        print("original bins",orig_bins)

        rmf_e=fits.open(self.input_rsp.rmf_path)[1]

        orig_rsp_fn="original_rsp_assembled.fits"
        fits.HDUList([fits.PrimaryHDU(),orig_e,rmf_e]).writeto(orig_rsp_fn,overwrite=True)

        ebins_compress_fn="ebins.txt"
        with open(ebins_compress_fn,"w") as f:
            for oi,(e1,e2) in enumerate(new_bins):
                i1 = list(orig_e.data['E_MIN']).index(e1)
                i2 = list(orig_e.data['E_MAX']).index(e2)
                print(oi,i1,i2,e1,e2)
                f.write("%i %i %i\n"%(i1,i2,i2-i1+1))

        new_rsp_fn="rsp_rebinned.fits"

        ht = pilton.heatool("rbnrmf")
        ht['infile']=orig_rsp_fn
        ht['binfile']=ebins_compress_fn
        ht['outfile']=new_rsp_fn
        ht['clobber']="yes"
        ht.run()


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
