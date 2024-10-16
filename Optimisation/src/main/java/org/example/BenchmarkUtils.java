package org.example;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.orekit.data.DataContext;
import org.orekit.data.DataProvider;
import org.orekit.data.DirectoryCrawler;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class BenchmarkUtils {

    public static void addData(final StringBuilder stringBuilder, final double timeFromStart, final double... data) {
        stringBuilder.append(timeFromStart);
        for (double element : data) {
            stringBuilder.append(",");
            stringBuilder.append(element);
        }
        stringBuilder.append("\n");
    }

   public static void addData(final StringBuilder stringBuilder, final double timeFromStart, final double[]... data
                             ) {
      stringBuilder.append(timeFromStart);
      for (double[] dataArray : data) {
         stringBuilder.append(",");
         stringBuilder.append(convertDoubleArrayToCSV(dataArray));
      }
      stringBuilder.append("\n");
   }

   public static void addData(final StringBuilder stringBuilder, final double[]... data) {
      for (double[] dataArray : data) {
         stringBuilder.append(convertDoubleArrayToCSV(dataArray));
         stringBuilder.append(",");
      }
      stringBuilder.append("\n");
   }

    public static void addData(final StringBuilder stringBuilder, final double... data) {
        for (double dataElement : data) {
            addDataElement(stringBuilder, dataElement);
        }
        stringBuilder.append("\n");
    }

   public static void addDataElement(final StringBuilder stringBuilder, final double data) {
      stringBuilder.append(data);
      stringBuilder.append(",");
   }

   public static String convertDoubleArrayToCSV(final double[] data) {

      return Stream.of(Arrays.toString(data).replaceAll("\\[(.*?)\\]", "$1"))
            .map(String::valueOf)
            .collect(Collectors.joining(","));
   }

   public static double[] convertRealMatrixToOneDimensionArray(final RealMatrix matrixToConvert) {

      final int rowDimension    = matrixToConvert.getRowDimension();
      final int columnDimension = matrixToConvert.getColumnDimension();

      final double[][] matrixData = matrixToConvert.getData();

      final double[] matrixDataOutput =
            new double[rowDimension * columnDimension];

      for (int row = 0; row < matrixToConvert.getRowDimension(); row++) {
          if (matrixToConvert.getColumnDimension() >= 0)
              System.arraycopy(matrixData[row], 0, matrixDataOutput, row * rowDimension,
                               matrixToConvert.getColumnDimension());
      }

      return matrixDataOutput;
   }

   public static void loadOrekitData() {

      //Loading Data
      final File rootFile   = new File(System.getProperty("user.home"));
      final File orekitData = new File(rootFile, "orekit-data");
      if (!orekitData.exists()) {
         System.out.format("Le fichier %s n'a pas été trouvé.%n", orekitData.getAbsolutePath());
      }
      final DataProvider dirCrawler = new DirectoryCrawler(orekitData);
      DataContext.getDefault().getDataProvidersManager().addProvider(dirCrawler);
   }

   public static void writeData(final File outputFile, final StringBuilder data) {

      try {
         final FileWriter fileWriter = new FileWriter(outputFile);
         fileWriter.write(data.toString());
         fileWriter.close();
      }
      catch (IOException e) {
         throw new RuntimeException(e);
      }
   }

   public static double computeRMS(final double[] data) {
      double sum = 0;
      for (final double element : data) {
         sum += element * element;
      }
      return FastMath.sqrt(sum / data.length);
   }
}
