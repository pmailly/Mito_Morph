package Mito_Utils;


import Mito_Stardist.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.ZProjector;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Object3D;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Objects3DPopulationColocalisation;
import mcib3d.geom.PairColocalisation;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mpicbg.ij.integral.RemoveOutliers;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;
import sc.fiji.analyzeSkeleton.AnalyzeSkeleton_;
import sc.fiji.analyzeSkeleton.SkeletonResult;


 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose dots_Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */
public class Mito_Processing {
    
    // DNA filter size
    private double minDNA = 0.004;
    private double maxDNA = 1;
    // Nucleus filter
    private double minNuc = 50;
    private double maxNuc = Double.MAX_VALUE;
    // Mito filter size
    private double minMito = 0.05;
    private double maxMito = Double.MAX_VALUE;
    
    // Max distance to track object
    public double maxDist = 5;
    private Calibration cal = new Calibration(); 
    public CLIJ2 clij2 = CLIJ2.getInstance();
    
    
    public String[] dotsDetectors = {"StarDist","DOG"};
    public String dotsDetector = "StarDist";
    
    // Stardist
    public Object syncObject = new Object();
    public final double stardistPercentileBottom = 0.2;
    public final double stardistPercentileTop = 99.8;
    public final double stardistProbThresh = 0.1;
    public final double stardistOverlayThresh = 0.35;
    public File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
    public String stardistModel = "";
    public String stardistOutput = "Label Image"; 
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    
    
     /**
     * check  installed modules
     * @return 
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /*
    Find starDist models in Fiji models folder
    */
    private String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        return(models);
    }
    
    
    /**
     * Find images in folder
     * @param imagesFolder
     * @param imageExt
     * @return 
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No Image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        return(images);
    }
       
     /**
     * Find channels name
     * @param imageName
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader, boolean bioformat) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                    if (!bioformat) {
                        channels[n] = channels[n].replace("_", "-");
                        channels[n] = "w"+(n+1)+channels[n];
                    }
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[0] = Integer.toString(n);
        }
        return(channels);         
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta, ImageProcessorReader reader) {
        cal = new Calibration();  
        // read image calibration
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        return(cal);
    }
    
    public Calibration getCalib()
    {
        System.out.println("x cal = " +cal.pixelWidth+", z cal=" + cal.pixelDepth);
        return cal;
    }
    
    
    /**
    * 
    * @param FileResults
    * @param resultsFileName
    * @param header
    * @return 
    */
    public BufferedWriter writeHeaders(String outFileResults, String header) throws IOException {
       FileWriter FileResults = new FileWriter(outFileResults, false);
       BufferedWriter outPutResults = new BufferedWriter(FileResults); 
       outPutResults.write(header);
       outPutResults.flush();
       return(outPutResults);
    }    
    
    /**
     * Dialog 
     */
    public int[] dialog(String[] channels) {
        String[] models = findStardistModels();
        String[] chNames = {"Nucleus", "DNA", "Mitos"};
        String dir = "";
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels selection", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String chName : channels) {
            gd.addChoice(chNames[index]+" : ", channels, channels[index]);
            index++;
        }
        gd.addMessage("DNA parameters", Font.getFont("Monospace"), Color.blue);
        gd.addMessage("Dots detection method", Font.getFont("Monospace"), Color.blue);
        gd.addChoice("Dots segmentation method :",dotsDetectors, dotsDetectors[0]);
        gd.addMessage("StarDist model", Font.getFont("Monospace"), Color.blue);
        if (models.length >= 2) {
            gd.addChoice("StarDist model :",models, models[0]);
        }
        else {
            gd.addMessage("No StarDist model found in Fiji !!", Font.getFont("Monospace"), Color.red);
            gd.addFileField("StarDist model :", stardistModel);
        }
        gd.addNumericField("Min DNA size (µm3) : ", minDNA, 3);
        gd.addNumericField("Max DNA size (µm3) : ", maxDNA, 3);
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Calibration xy (µm)  :", cal.pixelWidth, 3);
        if (cal.pixelDepth == 1)
            cal.pixelDepth = 0.5;
        gd.addNumericField("Calibration z (µm)  :", cal.pixelDepth, 3);
        gd.showDialog();
        int[] chChoices = new int[channels.length];
        for (int n = 0; n < chChoices.length; n++) {
            chChoices[n] = ArrayUtils.indexOf(channels, gd.getNextChoice());
        }
        dotsDetector = gd.getNextChoice();
        if (models.length >= 2) {
            stardistModel = modelsPath+File.separator+gd.getNextChoice();
        }
        else {
            stardistModel = gd.getNextString();
        }
        if (dotsDetector.equals("StarDist") && stardistModel.isEmpty()) {
            IJ.error("No model specify !!");
            return(null);
        }
        minDNA = gd.getNextNumber();
        maxDNA = gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelHeight = cal.pixelWidth;
        cal.pixelDepth = gd.getNextNumber();
        if (gd.wasCanceled())
                chChoices = null;
        return(chChoices);
    }
    
    
   /* Median filter 
     * Using CLIJ2
     * @param ClearCLBuffer
     * @param sizeXY
     * @param sizeZ
     */ 
    public ClearCLBuffer median_filter(ClearCLBuffer  imgCL, double sizeXY, double sizeZ) {
        ClearCLBuffer imgCLMed = clij2.create(imgCL);
        clij2.mean3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
        clij2.release(imgCL);
        return(imgCLMed);
    }
    
    /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ClearCLBuffer DOG(ClearCLBuffer imgCL, double size1, double size2) {
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        return(imgCLDOG);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     * @param fill 
     */
    public ClearCLBuffer threshold(ClearCLBuffer imgCL, String thMed) {
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        return(imgCLBin);
    }
    
    /**
     * Remove Outliers
     * 
     * @param img
     * @param radX
     * @param radY
     * @param factor
     * @return img
     */
    public void removeOutliers(ImagePlus img, int radX, int radY, float factor) {
        for (int i = 1; i <= img.getNSlices(); i++) {
            img.setSlice(i);
            ImageProcessor ip = img.getProcessor();
            RemoveOutliers removeOut = new RemoveOutliers(ip.convertToFloatProcessor());
            removeOut.removeOutliers(radX, radY, factor);
        }
    }

     /**
     * Nucleus segmentation
     * @param imgNuc
     * @param roi
     * @return 
     */
    public Objects3DPopulation find_nucleus(ImagePlus img, Roi roi) {
        ImagePlus img2 = img.duplicate();
        removeOutliers(img2, 50, 50, 1);
        ClearCLBuffer imgCL = clij2.push(img2);
        flush_close(img2);
        ClearCLBuffer imgCLBlur = clij2.create(imgCL);
        clij2.gaussianBlur2D(imgCL, imgCLBlur, 2, 2);
        clij2.release(imgCL);
        ClearCLBuffer imgCLBin = threshold(imgCLBlur, "Otsu");
        clij2.release(imgCLBlur);
        ClearCLBuffer imgFill = clij2.create(imgCLBin);
        clij2.binaryFillHoles(imgCLBin, imgFill);
        clij2.release(imgCLBin);
        ClearCLBuffer imgLabelled = clij2.create(imgFill);
        clij2.connectedComponentsLabelingBox(imgFill, imgLabelled);
        ImagePlus imgBin = clij2.pull(imgLabelled);
        clij2.release(imgLabelled);
        clearOutSide(imgBin, roi);
        imgBin.setCalibration(cal);
        Objects3DPopulation nucPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minNuc, maxNuc, true));
        flush_close(imgBin);
        return(nucPop);
    }
    
    /**
     * Mito segmentation
     * @param img
     * @return 
     */
    public Objects3DPopulation find_Mito(ImagePlus img, Roi roi, Objects3DPopulation nucPop) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgDOG = DOG(imgCL, 2, 3);
        clij2.release(imgCL);
        ClearCLBuffer imgCLBin = threshold(imgDOG, "Triangle");
        clij2.release(imgDOG);
        ClearCLBuffer imgLabelled = clij2.create(imgDOG);
        clij2.connectedComponentsLabelingBox(imgCLBin, imgLabelled);
        clij2.release(imgCLBin);
        ImagePlus imgBin = clij2.pull(imgLabelled);
        clij2.release(imgLabelled);
        imgBin.setCalibration(cal);
        clearOutSide(imgBin, roi);
        nucPop.draw(imgBin.getImageStack(), 0);
        Objects3DPopulation mitoPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minMito, maxMito, true));
        return(mitoPop);
    } 
    
    /**
     * DNA segmentation
     * @param img
     * @return 
     */
    public Objects3DPopulation find_DNA(ImagePlus img, Roi roi, Objects3DPopulation nucPop) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgDOG = DOG(imgCL, 2, 3);
        clij2.release(imgCL);
        ClearCLBuffer imgCLBin = threshold(imgDOG, "Triangle");
        clij2.release(imgDOG);
        ClearCLBuffer imgLabelled = clij2.create(imgCLBin);
        clij2.connectedComponentsLabelingBox(imgCLBin, imgLabelled);
        clij2.release(imgCLBin);
        ImagePlus imgBin = clij2.pull(imgLabelled);
        clij2.release(imgLabelled);
        imgBin.setCalibration(cal);
        clearOutSide(imgBin, roi);
        nucPop.draw(imgBin.getImageStack(), 0);
        Objects3DPopulation dnaPop = new Objects3DPopulation(getPopFromImage(imgBin).getObjectsWithinVolume(minDNA, maxDNA, true));
        flush_close(imgBin);
        return(dnaPop);
    }
    
    /**
     * Find gene population with Stardist
     */
    public Objects3DPopulation stardistDnaPop(ImagePlus img, Roi roi, Objects3DPopulation nucPop) throws IOException{
        // Go StarDist
        File starDistModelFile = new File(stardistModel);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(img);
        star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThresh, stardistOverlayThresh, stardistOutput);
        star.run();
        // label in 3D
        ImagePlus imgLab = star.getLabelImagePlus().duplicate();
        imgLab.setCalibration(cal);
        clearOutSide(imgLab, roi);
        nucPop.draw(imgLab.getImageStack(), 0);
        Objects3DPopulation dnaPop = new Objects3DPopulation(getPopFromImage(imgLab).getObjectsWithinVolume(minDNA, maxDNA, true));
        flush_close(imgLab);
        return(dnaPop);
        }
    
    
    public Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }
    
    
    // Threshold images and fill holes
    public void threshold(ImagePlus img, AutoThresholder.Method thMed, boolean fill, boolean calcul) {
        //  Threshold and binarize
       String cal = "";
       img.setZ(img.getNSlices()/2);
       img.updateAndDraw();
       IJ.setAutoThreshold(img, thMed.toString()+" dark");
       Prefs.blackBackground = false;
       if (calcul)
           cal = " calculate";
        IJ.run(img, "Convert to Mask", "method="+thMed.toString()+cal+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    // write object labels
    public void labelsObject (Object3D obj, ImagePlus img, int number, int color) {
        Font tagFont = new Font("SansSerif", Font.PLAIN, 32);
        int[] box = obj.getBoundingBox();
        int z = (int)obj.getCenterZ();
        int x = box[0] - 2;
        int y = box[2] - 2;
        img.setSlice(z+1);
        ImageProcessor ip = img.getProcessor();
        ip.setFont(tagFont);
        ip.setColor(color);
        ip.drawString(Integer.toString(number), x, y);
        img.updateAndDraw();    
    }
    
    // tag object number with random color
    public void tagsObject(ImageHandler imh, Object3D obj) {        
        int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
        obj.draw(imh, col);  
    }
    
   /**
     * ramdom color nucleus population
     */
    public ImagePlus randomColorPop (Objects3DPopulation cellsPop,  ImageHandler img, boolean label) {
        //create image objects population
        img.set332RGBLut();
        img.setCalibration(img.getCalibration());
        for (int i = 0; i < cellsPop.getNbObjects(); i++) {
            Object3D obj = cellsPop.getObject(i);
            int col = ThreadLocalRandom.current().nextInt(2, 255 + 1);
            obj.draw(img, col);
            if (label)
               labelsObject(obj, img.getImagePlus(), (i+1), col); 
        } 
        return(img.getImagePlus());
    } 

 
    public ImagePlus doZProjection(ImagePlus img) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(ZProjector.MAX_METHOD);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
// Flush and close images
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
    
    
    /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        roi.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(roi);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(roi);
        }
        img.deleteRoi();
        img.updateAndDraw();
    }
    
    
    /**
     * Analayze skeleton
     * @param img
     * @param output
     * @return {#branch, branchLenght, #endPoint, #junction}
     */

    public double[] analyzeSkeleton (ImagePlus img, Roi roi, Objects3DPopulation pop, String outFileName) {
        ImageHandler imh = ImageHandler.wrap(img).createSameDimensions();
        pop.draw(imh, 255);
        ImagePlus imgBin = imh.getImagePlus();
        imh.closeImagePlus();
        imgBin.setCalibration(cal);
        IJ.run(imgBin,"8-bit","");
        IJ.run(imgBin, "Skeletonize (2D/3D)", "");
	String imgTitle = img.getTitle();
        AnalyzeSkeleton_ analyzeSkeleton = new AnalyzeSkeleton_();
        AnalyzeSkeleton_.calculateShortestPath = true;
        analyzeSkeleton.setup("",imgBin);
        SkeletonResult skeletonResults = analyzeSkeleton.run(AnalyzeSkeleton_.NONE, false, true, null, true, false);
        ImageStack imgStackLab = analyzeSkeleton.getLabeledSkeletons();
        //  compute parameters for each skeleton
        IJ.showStatus("Computing parameters for each skeleton ...");
        ImagePlus imgLab = new ImagePlus(imgTitle+"_LabelledSkel.tif", imgStackLab);
        ImagePlus imgLabProj = doZProjection(imgLab);
        IJ.run(imgLabProj, "3-3-2 RGB", "");
        imgLabProj.setCalibration(img.getCalibration());
        flush_close(imgLab);
        int[] branchNumbers = skeletonResults.getBranches();
        double[] branchLengths = skeletonResults.getAverageBranchLength();
        int[] nbEndPoints = skeletonResults.getEndPoints();
        int[] junctions = skeletonResults.getJunctions();
        int branches = 0;
        double branchLength = 0;
        int endPoint = 0;
        int junction = 0;
        for (int i = 0; i < skeletonResults.getGraph().length; i++) {
            branches += branchNumbers[i];
            branchLength += branchLengths[i];
            endPoint += nbEndPoints[i];
            junction += junctions[i];
        }
        double[] skeletonParams = {branches, branchLength, endPoint, junction};
        FileSaver imgSave = new FileSaver(imgLabProj);
        imgSave.saveAsTiff(outFileName+"_Roi-"+roi.getName()+"_LabelledSkel.tif");
        flush_close(imgLabProj); 
        flush_close(imgBin);
        return(skeletonParams);
    }
    
    /**
     * find dna inside Mito
     * @param dnaPop
     * @param mitoPop
     * @return 
     */
    public Objects3DPopulation findDnaInMito(Objects3DPopulation dnaPop, Objects3DPopulation mitoPop) {
        Objects3DPopulation dnaMitoPop = new Objects3DPopulation();
        Objects3DPopulationColocalisation coloc =  new Objects3DPopulationColocalisation(dnaPop, mitoPop);
        ArrayList<PairColocalisation> pairColoc = coloc.getAllColocalisationPairs();
        for (PairColocalisation p : pairColoc) {
            int volColoc = p.getVolumeColoc();
            if (volColoc > 0)
                dnaMitoPop.addObject(p.getObject3D1());
        }
        dnaMitoPop.setCalibration(cal.pixelWidth, cal.pixelDepth, cal.getUnit());
        return(dnaMitoPop);
    }

    
    /**
    * Local compute parameters add to results file
    * @param nuc nucleus number 
    * @param mitoPop mito population
    * @param mitoParams branch number, branch lenght, end points, junctions
     * @param imgName
    * @param results buffer
    * @throws java.io.IOException
    **/
    public void computeParameters(int nuc, Objects3DPopulation mitoPop, ImagePlus imgMito, Objects3DPopulation dnaPop, ImagePlus imgDna, 
            Objects3DPopulation dnaInMitoPop, double[] mitoParams, String roi, 
            String imgName, BufferedWriter results) throws IOException {
        IJ.showStatus("Computing parameters ....");
        // mito volume
        double mitoVol = 0, dnaVol = 0, dnaInMitoVol = 0;
        double mitoInt = 0, dnaInt = 0, dnaInMitoInt = 0;
        int mitos = mitoPop.getNbObjects();
        for (int i = 0; i < mitos; i++) {
            mitoVol += mitoPop.getObject(i).getVolumeUnit();
            mitoInt += mitoPop.getObject(i).getIntegratedDensity(ImageHandler.wrap(imgMito));
        }
        int dnas = dnaPop.getNbObjects();
        for (int i = 0; i < dnas; i++) {
            dnaVol += dnaPop.getObject(i).getVolumeUnit();
            dnaInt += dnaPop.getObject(i).getIntegratedDensity(ImageHandler.wrap(imgDna));
        }
        int dnaInMito = dnaInMitoPop.getNbObjects();
        for (int i = 0; i < dnaInMito; i++) {
            dnaInMitoVol += dnaInMitoPop.getObject(i).getVolumeUnit();
            dnaInMitoInt += dnaInMitoPop.getObject(i).getIntegratedDensity(ImageHandler.wrap(imgDna));
        }
        
        results.write(imgName+"\t"+roi+"\t"+nuc+"\t"+mitos+"\t"+mitoVol+"\t"+mitoInt+"\t"+mitoParams[0]+"\t"+mitoParams[1]+"\t"+mitoParams[2]+"\t"+
                mitoParams[3]+"\t"+dnas+"\t"+dnaVol+"\t"+dnaInt+"\t"+dnaInMito+"\t"+dnaInMitoVol+"\t"+dnaInMitoInt
                +"\t"+(dnas - dnaInMito)+"\t"+(dnaVol - dnaInMitoVol)+"\t"+(dnaInt -dnaInMitoInt)+"\n");
        results.flush();
    }
}
