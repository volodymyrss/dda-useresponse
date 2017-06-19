import glob

import ddosa
import findarf

def test_find():
    fa=findarf.FindARF(assume=ddosa.ScWData(input_scwid="066500220010.001"))
    fa.get()

    assert fa.arf_path == "/home/isdc/savchenk/osa11_deployment/deployment/ic/ic/ibis/rsp/isgr_arf_rsp_0035.fits"
