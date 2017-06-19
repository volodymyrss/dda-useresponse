import glob

import ddosa
import useresponse

def test_find():
    ur=useresponse.FindResponse(assume=ddosa.ScWData(input_scwid="066500220010.001"))
    ur.get()

    assert ur.arf_path == "/home/isdc/savchenk/osa11_deployment/deployment/ic/ic/ibis/rsp/isgr_arf_rsp_0037.fits"
