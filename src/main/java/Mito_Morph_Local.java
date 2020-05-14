/*
 * Analyze mitochondrial morphological network
 * 
 * Author Philippe Mailly
 */


/* 
* Images on local machine
*/


import static Mito_Utils.Mito_Processing.analyzeSkeleton;
import static Mito_Utils.Mito_Processing.computeMitoParameters;
import static Mito_Utils.Mito_Processing.find_nucleus2;
import static Mito_Utils.Mito_Processing.getPopFromImage;
import ij.*;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.AutoThresholder;
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
import mcib3d.image3d.ImageInt;
import static Mito_Utils.Mito_Processing.flush_close;
import static Mito_Utils.Mito_Processing.median_filter;
import static Mito_Utils.Mito_Processing.randomColorPop;
import static Mito_Utils.Mito_Processing.threshold;
import static Mito_Utils.JDialogOmeroConnect.dialogCancel;
import static Mito_Utils.JDialogOmeroConnect.imagesFolder;
import static Mito_Utils.Mito_Processing.clearOutSide;
import static Mito_Utils.Mito_Processing.writeHeaders;
import ij.gui.Roi;
import ij.plugin.RGBStackMerge;
import ij.plugin.frame.RoiManager;
import java.util.Arrays;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import org.apache.commons.io.FilenameUtils;


public class Mito_Morph_Local implements PlugIn {

    private String imageDir = "";
    public static String outDirResults = "";
    public static final Calibration cal = new Calibration();


// min volume in microns^3 for dots
    private final double minMito = 0.05;
// max volume in microns^3 for dots
    private final double maxMito = Double.MAX_VALUE;
    public BufferedWriter outPutGlobalResults;    
    
    

    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        String imageExt = "czi";
        final boolean canceled = false;
        
        try {
            if (dialogCancel){
                IJ.showStatus(" Pluging canceled");
                return;
            }
            if (imagesFolder == null) {
                return;
            }
            File inDir = new File(imagesFolder);
            String[] imageFile = inDir.list();
            if (imageFile == null) {
                System.out.println("No Image found in "+imagesFolder);
                return;
            }
            // create output folder
            outDirResults = imagesFolder + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }

            // create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            Arrays.sort(imageFile);
            int imageNum = 0; 
            for (int i = 0; i < imageFile.length; i++) {
                String fileExt = FilenameUtils.getExtension(imageFile[i]);
                if (fileExt.equals(imageExt)) {
                    int series = 0;
                    String imageName = inDir+ File.separator+imageFile[i];
                    String rootName = FilenameUtils.getBaseName(imageFile[i]);
                    imageNum++;
                    reader.setId(imageName);
                    // Check calibration
                    if (imageNum == 1) {
                        cal.pixelWidth = meta.getPixelsPhysicalSizeX(series).value().doubleValue();
                        cal.pixelHeight = cal.pixelWidth;
                        if (meta.getPixelsPhysicalSizeZ(series) != null)
                            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(series).value().doubleValue();
                        else
                            cal.pixelDepth = 1;
                        cal.setUnit("microns");
                        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);

                        /*
                        * Write headers results for results file
                        */
                        // Global file for mito results
                        String resultsName = "GlobalResults.xls";
                        String header = "ImageName\tRoi\tNucleus number\tMito number\tMito Volume\tMito branch number\tMito branch length\t"
                                + "Mito end points\tMito junction number\n";
                        outPutGlobalResults = writeHeaders(outDirResults+resultsName, header); 
                    }

                    series = reader.getSeriesCount();  
                    for (int s = 0; s < series; s++) {
                        reader.setSeries(s);
                        String seriesName = meta.getImageName(s);

                        ImporterOptions options = new ImporterOptions();
                        options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                        options.setId(imageName);
                        options.setSplitChannels(true);
                        options.setQuiet(true);
                        options.setSeriesOn(s, true);
                        
                        // find mito only inside roi if exist roi file
                        String roiFile = "";
                        if (new File(inDir + File.separator+rootName + ".roi").exists())
                            roiFile = inDir + File.separator+rootName + ".roi";
                        else if (new File(inDir + File.separator+rootName + ".zip").exists())
                            roiFile = inDir + File.separator+rootName + ".zip";
                        Roi[] rois = null;
                        if (!"".equals(roiFile)) {
                            // find rois
                            System.out.println("Find roi " + new File(roiFile).getName());
                            RoiManager rm = new RoiManager(false);
                            rm.runCommand("Open", roiFile);
                            rois = rm.getRoisAsArray();
                        }
                        else {
                            //take the whole image
                            IJ.showMessage("No roi found, taking the whole image");
                            rois[0] = new Roi(0, 0, reader.getSizeX(), reader.getSizeY());
                        }
                        /*
                        * Open channels
                        */
                        
                        // DAPI
                        int dapiCh = 1;                        
                        options.setCBegin(s, dapiCh);
                        options.setCEnd(s, dapiCh);
                        System.out.println("-- Series : "+ seriesName);
                        System.out.println("Opening Nucleus channel");
                        ImagePlus imgNucOrg= BF.openImagePlus(options)[0]; 

                        // mito channel
                        int mitoCh = 0;
                        System.out.println("Opening mito channel");
                        options.setCBegin(s, mitoCh);
                        options.setCEnd(s, mitoCh);
                        ImagePlus imgMitoOrg = BF.openImagePlus(options)[0];
                        
                        // for each roi
                        
                        for( int r = 0; r < rois.length; r++) {
                            Roi roi = rois[r];
                            
                            // Find nucleus
                            Objects3DPopulation nucPop = new Objects3DPopulation();
                            imgNucOrg.setRoi(roi);
                            ImagePlus imgNuc = imgNucOrg.duplicate();
                            nucPop = find_nucleus2(imgNuc, roi);
                            int totalNucPop = nucPop.getNbObjects();
                            System.out.println("Roi " + (r+1)+" Detected nucleus = "+totalNucPop);

                            // Find mitos
                            imgMitoOrg.setRoi(roi);
                            ImagePlus imgMito = imgMitoOrg.duplicate();
                            median_filter(imgMito, 1.5);
                            IJ.run(imgMito, "Laplacian of Gaussian", "sigma=4 scale_normalised negate stack");
                            threshold(imgMito, AutoThresholder.Method.RenyiEntropy, false, false);
                            clearOutSide(imgMito, roi);

                            Objects3DPopulation mitoPop = getPopFromImage(imgMito, cal);
                            //objectsSizeFilter(minMito, maxMito, mitoPop, imgMitoOrg, false); 
                            System.out.println("Mito pop = "+ mitoPop.getNbObjects());

                            // Find mito network morphology
                            // Skeletonize
                            double[] skeletonParams = analyzeSkeleton(imgMito, outDirResults+rootName);


                            // Compute global Mito parameters                        
                            // nb of mito, mean mito volume, skeleton parameters
                            IJ.showStatus("Writing parameters ...");
                            computeMitoParameters(nucPop.getNbObjects(), mitoPop, skeletonParams, (r+1), rootName+seriesName+"_Mito", outPutGlobalResults);

                            // Save objects image
                            ImageHandler imhMitoObjects = ImageHandler.wrap(imgMito).createSameDimensions();
                            ImageHandler imhNucObjects = imhMitoObjects.duplicate();
                            mitoPop.draw(imhMitoObjects, 255);
                            nucPop.draw(imhNucObjects, 255);
                            ImagePlus[] imgColors = {imhMitoObjects.getImagePlus(), null, imhNucObjects.getImagePlus()};
                            ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
                            imgObjects.setCalibration(cal);
                            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
                            FileSaver ImgObjectsFile = new FileSaver(imgObjects);
                            ImgObjectsFile.saveAsTiff(outDirResults + rootName + "_" + seriesName + "Roi_"+(r+1)+"_Objects.tif");
                            flush_close(imgObjects);
                            flush_close(imhMitoObjects.getImagePlus());
                            flush_close(imhNucObjects.getImagePlus());
                            flush_close(imgNuc);
                            flush_close(imgMito);
                            options.setSeriesOn(s, false);
                        }
                        flush_close(imgNucOrg);
                        flush_close(imgMitoOrg);
                    }
                }
            }
            outPutGlobalResults.close();
            IJ.showStatus("Process done");
        }   catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Mito_Morph.class.getName()).log(Level.SEVERE, null, ex);
        } 
    }
}