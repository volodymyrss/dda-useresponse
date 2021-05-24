import glob
import os

import astropy.io.fits as fits
import numpy as np

import ast
import ddosa
import dataanalysis.core as da
import dataanalysis.hashtools as ht
import pilton

import dataanalysis.importing as importing
import pilton
import findic

try:
    import findic
except:
    findic,findic_name=importing.load_by_name("git://findic")

class FindICARF(findic.FindICIndexEntry):
    ds="ISGR-ARF.-RSP"


class FindICRMF(findic.FindICIndexEntry):
    ds="ISGR-RMF.-RSP"

class ICRMF(ddosa.DataAnalysis):
    pass

class FindResponse(ddosa.DataAnalysis):
#    input_findicarf = FindICARF
    input_scw = ddosa.ScWData
    input_findicrmf = FindICRMF
    input_icroot = ddosa.ICRoot

    version="v2"

    #run_for_hashe=True

    def main(self):
        self.ic_entry = self.input_findicrmf.find_entry(self.input_scw, self.input_icroot)
        self.rmf_path = self.ic_entry['member_location']

        print(self.ic_entry)

        return ICRMF(input_icroot=self.input_icroot, use_rmf_path=self.rmf_path, version="v1-"+self.ic_entry['idx_hash'])

class FindICEBDS(findic.FindICIndexEntry):
    ds="ISGR-EBDS-MOD"

    input_scw=da.NoAnalysis


class CompressEBins(ddosa.DataAnalysis):
    input_ic_ebds=FindICEBDS

    factor=1

    def get_version(self):
        v=super().get_version()
        v+="_f%i"%self.factor
        return v

    def main(self):
        self.ic_ebds_member_location = self.input_ic_ebds.get_member_location()
        e = fits.open(self.ic_ebds_member_location)[1].data
        ic_bins = zip(e['E_MIN'], e['E_MAX'])

        print("raw bins",ic_bins)

        o_e1=e['E_MIN'][::int(self.factor)]
        o_e2=e['E_MAX'][self.factor-1::int(self.factor)]

        o_ebins=zip(o_e1,o_e2)
        o_e1,o_e2=zip(*o_ebins)

        print("compressed bins",o_ebins)

        f=fits.open(self.ic_ebds_member_location)
        f[1].data=f[1].data[:len(o_e1)]

        print(len(f[1].data['E_MIN']),len(o_e1))

        self.ebds_mod_fn = "compressed_ebins.fits"

        f[1].data['E_MIN'] = o_e1
        f[1].data['E_MAX'] = o_e2
        f.writeto(self.ebds_mod_fn,overwrite=True)


class RebinResponse(ddosa.DataAnalysis):
    input_rsp = FindResponse
    #input_ic_ebds = FindICEBDS
    input_ebins = CompressEBins

    version="v1"

    cached=True

    @property
    def rmf_path(self):
        return self.rmf.get_path()

    cached=False

    def main(self):
        new_e = fits.open(self.input_ebins.ebds_mod_fn)[1]
        new_bins = zip(new_e.data['E_MIN'], new_e.data['E_MAX'])
        print("new bins:",new_bins)

        orig_e = fits.open(self.input_ebins.ic_ebds_member_location)[1]
        orig_bins = zip(orig_e.data['E_MIN'], orig_e.data['E_MAX'])
        print("original bins",orig_bins)

        print("original rmf",self.input_rsp.rmf_path)

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

        try:
            ht = pilton.heatool("rbnrmf")
            ht['infile'] = orig_rsp_fn
            ht['binfile'] = ebins_compress_fn
            ht['outfile'] = new_rsp_fn
            ht['clobber'] = "yes"
            ht.run()
        except Exception as e:
            print("problem running rbnrmf", e)
            raise

        self.rmf=da.DataFile(new_rsp_fn)

    @property
    def rmf_path(self):
        return self.rmf.get_path()

class RebinResponseProcessingSummary(ddosa.DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=RebinResponse(assume=[ddosa.ScWData(input_scwid="any",use_abstract=True),ddosa.Revolution(input_revid="any",use_abstract=True)]) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
        print("generalized hash:",ahash)
        rh=ht.shhash(ahash)
        print("reduced hash",rh)
        handle=da.DataHandle('processing_definition:'+rh[:8])
        self.factory.note_factorization(dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=ahash,
            reduced_hash=rh,
            handle=handle.handle,
        ))
        return [handle]

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
