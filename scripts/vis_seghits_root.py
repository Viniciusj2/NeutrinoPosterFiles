import ROOT as rt
from math import sqrt

def vis_seghits_root( seghit_v, cryostat_radius_cm=37.50, cryostat_height_cm=133.70, entry=None ):

    cvis = rt.TCanvas("cvisseghits", "Visualizae Seghits", 1500,600)
    cvis.Divide(2,1)

    # just setting up the view for the event
    hxy = rt.TH2D("hxy","; x (cm); y (cm)",10, -3*cryostat_radius_cm, 3*cryostat_radius_cm, 10, -3*cryostat_radius_cm, 3*cryostat_radius_cm )
    hrz = rt.TH2D("hrz","; r (cm); z (cm)",10, 0, 3*cryostat_radius_cm, 10, -3*cryostat_height_cm/2.0, 3*cryostat_height_cm/2.0 )

    if entry is not None:
        hxy.SetTitle(entry)
        hrz.SetTitle(entry)

    # XY
    cvis.cd(1)
    hxy.Draw()
    tcirc = rt.TEllipse( 0.0, 0.0, cryostat_radius_cm, cryostat_radius_cm )
    tcirc.SetFillStyle(0)
    tcirc.SetLineColor( rt.kBlack )
    tcirc.Draw()
    
    gxy= rt.TGraph( seghit_v.size() )
    for ii in range(seghit_v.size()):
        hit = seghit_v.at(ii)
        gxy.SetPoint(ii, hit.GetStart()[0]*0.1, hit.GetStart()[1]*0.1 )
    gxy.SetMarkerColor( rt.kRed )
    gxy.SetMarkerStyle(20)
    gxy.Draw("P")
    cvis.Update()

    # RZ
    cvis.cd(2)
    hrz.Draw()
    tbox = rt.TBox( 0, -0.5*cryostat_height_cm, cryostat_radius_cm, 0.5*cryostat_height_cm )
    tbox.SetFillStyle(0)
    tbox.SetLineColor( rt.kBlack )
    tbox.Draw()
    
    grz= rt.TGraph( seghit_v.size() )
    for ii in range(seghit_v.size()):
        hit = seghit_v.at(ii)
        x = hit.GetStart()[0]*0.1
        y = hit.GetStart()[1]*0.1
        r = sqrt(x*x+y*y)
        z = hit.GetStart()[2]*0.1
        grz.SetPoint(ii, r, z )
    grz.SetMarkerColor( rt.kRed )
    grz.SetMarkerStyle(20)
    grz.Draw("P")
    cvis.Update()

    print("[ENTER] to continue")
    input()
    return [cvis, hxy, hrz, gxy, grz, tcirc, tbox]
    
    
    

