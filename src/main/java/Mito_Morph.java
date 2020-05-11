/*
 * Analyze mitochondrial morphological network
 * 
 * Author Philippe Mailly
 */

import ij.*;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import Mito_Utils.JDialogOmeroConnect;
import static Mito_Utils.JDialogOmeroConnect.localImages;
import java.awt.Frame;


public class Mito_Morph implements PlugIn {

    private final boolean dialogCancel = false;
    private final boolean canceled = false;
    private String imageDir = "";
    public static String outDirResults = "";
    public static final Calibration cal = new Calibration();


// min volume in microns^3 for dots
    private final double minMito = 0.05;
// max volume in microns^3 for dots
    private final double maxMito = Double.MAX_VALUE;
   
// Default Z step
    public static double zStep = 0.193;

    public BufferedWriter outPutGlobalResults;    
    
    public static String imageExt = "czi";

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        if (canceled) {
            IJ.showMessage(" Pluging canceled");
            return;
        }
        JDialogOmeroConnect dialog = new JDialogOmeroConnect(new Frame(), true);
        dialog.show();
        if (dialogCancel){
            IJ.showStatus(" Pluging canceled");
            return;
        }

        /* 
        * Images on local machine
        */

        if (localImages) {
            new Mito_Morph_Local().run("");
        }
        /*
        Images on OMERO server
        */

        else {
            new Mito_Morph_Omero().run("");     
        }

        IJ.showStatus("Process done");
    }
}