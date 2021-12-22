/*
 * Analyze mitochondrial morphological network
 * 
 * Author Philippe Mailly
 */


/* 
* Images on local machine
*/


import Mito_Utils.Mito_Processing;

import ij.*;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.image3d.ImageHandler;

import ij.gui.Roi;
import ij.plugin.RGBStackMerge;
import ij.plugin.frame.RoiManager;
import java.util.ArrayList;
import java.util.Collections;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import java.awt.Rectangle;
import loci.common.Region;
import org.apache.commons.io.FilenameUtils;


public class Mito_Morph implements PlugIn {

    private String imageDir = "";
    public static String outDirResults = "";
    private static Calibration cal = new Calibration();
    private File inDir;
    private String[] chsName;
    public BufferedWriter outPutGlobalResults;    
    
    Mito_Processing proc = new Mito_Processing();

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        String imageExt = "czi";
        final boolean canceled = false;
        
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            if (!proc.checkInstalledModules()) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Images folder");
            if (imageDir == null) {
                return;
            }
            inDir = new File(imageDir);
            ArrayList<String> imageFiles = proc.findImages(imageDir, "czi");
            if (imageFiles == null) {
                return;
            }
            // create output folder
            outDirResults = imageDir + "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            

            // Reset foreground and background
            IJ.run("Colors...", "foreground=white background=black");
            
            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Collections.sort(imageFiles);
            // Find channel names , calibration
            reader.setId(imageFiles.get(0));
            cal = proc.findImageCalib(meta, reader);
            chsName = proc.findChannels(imageFiles.get(0), meta, reader, true);
            int[] channelIndex = proc.dialog(chsName);
            cal = proc.getCalib();
            if (channelIndex == null)
                return;
            /*
            * Write headers results for results file
            */
            // Global file for mito results
            String resultsName = "GlobalResults-"+proc.dotsDetector+".xls";
            String header = "ImageName\tRoi\tNucleus number\tMito number\tMito volume\tMito intensity\tMito branch number\tMito branch length\t"
                    + "Mito end points\tMito junction number\tTotal DNA number\tTotal DNA volume\tTotal DNA mean intensity\tDNA in mito number\tDNA in mito volume"
                    + "\tDNA in mito mean intensity\tDNA out Mito number\tDNA out Mito volume\tDNA out mito mean intensity\n";
            outPutGlobalResults = proc.writeHeaders(outDirResults+resultsName, header); 
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                reader.setId(f);
                int series = reader.getSeries();
                reader.setSeries(series);
                series = reader.getSeriesCount();  
                for (int s = 0; s < series; s++) {
                    reader.setSeries(s);
                    String seriesName = meta.getImageName(s);
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setQuiet(true);
                    options.setSeriesOn(s, true);
                    options.setCrop(true);
                    // find rois if exist roi file
                    String roiFile = "";
                    if (new File(inDir + File.separator+rootName + ".roi").exists())
                        roiFile = inDir + File.separator+rootName + ".roi";
                    else if (new File(inDir + File.separator+rootName + ".zip").exists())
                        roiFile = inDir + File.separator+rootName + ".zip";
                    else {
                        IJ.showStatus("No roi file found skip image !!!");
                        return;
                    }
                    ArrayList<Roi> rois = new ArrayList<>();
                    // find rois
                    System.out.println("Find roi " + new File(roiFile).getName());
                    RoiManager rm = new RoiManager(false);
                    rm.runCommand("Open", roiFile);
                    for (Roi roi : rm.getRoisAsArray()) 
                        rois.add(roi);
                    
                    

                    // for each roi
                    for( int r = 0; r < rois.size(); r++) {
                        Roi roi = rois.get(r);
                        String roiName = roi.getName();
                        Rectangle rect = roi.getBounds();
                        Region reg = new Region(rect.x, rect.y, rect.width, rect.height);
                        options.setCropRegion(s, reg);
                        
                        /*
                        * Open channels
                        */
                        // DAPI channel
                        int dapiCh = channelIndex[0];                        
                        options.setCBegin(s, dapiCh);
                        options.setCEnd(s, dapiCh);System.out.println("-- Series : "+ seriesName);
                        System.out.println("Opening Nucleus channel");
                        ImagePlus imgNucOrg= BF.openImagePlus(options)[0];  
                        
                        // Find nucleus
                        Objects3DPopulation nucPop = proc.find_nucleus(imgNucOrg, roi);
                        int totalNucPop = nucPop.getNbObjects();
                        System.out.println(roiName +" roi Detected nucleus = "+totalNucPop);
                        
                        // Dna channel
                        int dnaCh = channelIndex[1];
                        options.setCBegin(s, dnaCh);
                        options.setCEnd(s, dnaCh);System.out.println("Opening dna channel");
                        ImagePlus imgDnaOrg = BF.openImagePlus(options)[0];
                         // Find Dna
                        Objects3DPopulation dnaPop = new Objects3DPopulation();
                        if (!proc.dotsDetector.equals("StarDist"))
                            dnaPop = proc.find_DNA(imgDnaOrg, roi, nucPop);
                        else
                            dnaPop = proc.stardistDnaPop(imgDnaOrg, roi, nucPop);
                        System.out.println("DNA pop = "+ dnaPop.getNbObjects());                        
                        
                        // Mito channel
                        int mitoCh = channelIndex[2];
                        options.setCBegin(s, mitoCh);
                        options.setCEnd(s, mitoCh);System.out.println("Opening mito channel");
                        ImagePlus imgMitoOrg = BF.openImagePlus(options)[0];// Find Mitos
                        Objects3DPopulation mitoPop = proc.find_Mito(imgMitoOrg, roi, nucPop);
                        System.out.println("Mito pop = "+ mitoPop.getNbObjects());
                        
                        // Find dna inside mito
                        Objects3DPopulation dnaInMitoPop = proc.findDnaInMito(dnaPop, mitoPop);
                        System.out.println("Dna inside mito = "+dnaInMitoPop.getNbObjects());
                        
                        // Find mito network morphology
                        // Skeletonize
                        double[] skeletonParams = proc.analyzeSkeleton(imgMitoOrg, roi, mitoPop, outDirResults+rootName);
                        // Compute global parameters                        
                        IJ.showStatus("Writing parameters ...");
                        proc.computeParameters(nucPop.getNbObjects(), mitoPop, imgMitoOrg, dnaPop, imgDnaOrg, dnaInMitoPop, skeletonParams, roiName, rootName, outPutGlobalResults);
                        proc.flush_close(imgMitoOrg);
                        proc.flush_close(imgDnaOrg);
                        
                        
                        // Save objects image
                        ImageHandler imhMitoObjects = ImageHandler.wrap(imgNucOrg).createSameDimensions();
                        ImageHandler imhNucObjects = imhMitoObjects.duplicate();
                        ImageHandler imhDnaObjects = imhMitoObjects.duplicate();
                        mitoPop.draw(imhMitoObjects, 255);
                        nucPop.draw(imhNucObjects, 255);
                        dnaPop.draw(imhDnaObjects, 255);
                        ImagePlus[] imgColors = {imhMitoObjects.getImagePlus(), imhDnaObjects.getImagePlus(), imhNucObjects.getImagePlus()};
                        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                        imgObjects.setCalibration(cal);
                        IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                        ImgObjectsFile.saveAsTiff(outDirResults + rootName + "_Roi-"+roiName+"_Objects-"+proc.dotsDetector+".tif");
                        proc.flush_close(imgObjects);
                        proc.flush_close(imhMitoObjects.getImagePlus());
                        proc.flush_close(imhNucObjects.getImagePlus());
                        proc.flush_close(imgNucOrg);
                    }
                    options.setSeriesOn(s, false);
                }
            }
            outPutGlobalResults.close();
            IJ.showStatus("Process done");
        }   catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Mito_Morph.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}