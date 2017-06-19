import glob
import os

import astropy.io.fits as fits
import numpy as np

import ddosa

class FindResponse(ddosa.DataAnalysis):
    input_scw=ddosa.ScWData
    input_ic=ddosa.ICRoot

    def main(self):
        t1,t2=self.input_scw.get_t()

        idx=fits.open(self.input_ic.icroot+"/idx/ic/ISGR-ARF.-RSP-IDX.fits")[1].data

        vstart=idx['VSTART']
        vstop=idx['VSTOP']

        m_v=(vstart<t1) & (vstop>t2)

        print t1,t2
        print idx[m_v]

        i_s=np.argsort(idx['VSTART'][m_v])
        print idx[m_v][i_s]
        i=i_s[-1]
        print idx[m_v][i]
        print vstart[m_v][i]-t1,vstop[m_v][i]-t2

        arf_path=idx['MEMBER_LOCATION'][m_v][i]
        arf_path_abs=os.path.abspath(self.input_ic.icroot+"/idx/ic/"+arf_path)
        print arf_path_abs

        print fits.open(arf_path_abs)[1].header

        self.arf_path=arf_path_abs
