/*
 * The MIT License
 *
 * Copyright (c) 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
package utils;

import casekit.nmr.model.Assignment;
import casekit.nmr.model.ExtendedConnectionMatrix;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.reflect.TypeToken;
import model.Query;
import model.SSC;
import org.bson.Document;

import java.io.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Consumer;

/**
 * @author Michael Wenk [https://github.com/michaelwenk]
 */
public class Converter {

    private final static Gson GSON = new GsonBuilder().setLenient()
                                                      .create(); //.setPrettyPrinting()

    public static String documentToJson(final Document document) {
        return document.toJson();
    }

    public static String sscToJson(final SSC ssc, final long sscCount) {
        //        return Converter.documentToJson(Converter.sscToDocument(ssc));
        return Converter.documentToJson(new Document().append(String.valueOf(sscCount), Converter.sscToDocument(ssc)));
    }

    public static Document sscToDocument(final SSC ssc) {
        final Document document = new Document();
        document.append("structure", Document.parse(GSON.toJson(
                GSON.toJsonTree(new ExtendedConnectionMatrix(ssc.getStructure()), ExtendedConnectionMatrix.class))));
        document.append("spectrum", Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getSpectrum(), Spectrum.class))));
        document.append("assignment",
                        Document.parse(GSON.toJson(GSON.toJsonTree(ssc.getAssignment(), Assignment.class))));
        document.append("hoseCodes", ssc.getHoseCodes());
        document.append("unsaturatedAtomIndices", ssc.getUnsaturatedAtomIndices());

        return document;
    }


    public static SSC jsonToSsc(final String json) {
        final JsonObject jsonObject = JsonParser.parseString(json)
                                                .getAsJsonObject();
        return new SSC(GSON.fromJson(jsonObject.get("structure"), ExtendedConnectionMatrix.class)
                           .toAtomContainer(), GSON.fromJson(jsonObject.get("spectrum"), Spectrum.class),
                       GSON.fromJson(jsonObject.get("assignment"), Assignment.class),
                       GSON.fromJson(jsonObject.get("hoseCodes"), new TypeToken<List<String>>() {
                       }.getType()),
                       GSON.fromJson(jsonObject.get("unsaturatedAtomIndices"), new TypeToken<List<Integer>>() {
                       }.getType()));
    }


    public static List<SSC> jsonToSscList(final String json, final int nThreads) throws InterruptedException {
        return Converter.jsonToSscList(new BufferedReader(new StringReader(json)), nThreads);
    }

    public static List<SSC> jsonFileToSscList(final String pathToJsonFile,
                                              final int nThreads) throws InterruptedException, FileNotFoundException {
        return Converter.jsonToSscList(new BufferedReader(new FileReader(pathToJsonFile)), nThreads);
    }

    private static List<SSC> jsonToSscList(final BufferedReader br, final int nThreads) throws InterruptedException {
        System.out.println("started converting... jsonToSscList");

        final ConcurrentLinkedQueue<SSC> convertedSSCs = new ConcurrentLinkedQueue<>();
        final List<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        br.lines()
          .forEach(line -> {
              if ((line.trim()
                       .length()
                      > 1)
                      || (!line.trim()
                               .startsWith("{")
                      && !line.trim()
                              .endsWith("}"))) {
                  final StringBuilder sscInJSON = new StringBuilder();
                  if (line.endsWith(",")) {
                      sscInJSON.append(line, 0, line.length()
                              - 1);
                  } else {
                      sscInJSON.append(line);
                  }
                  callables.add(() -> Converter.jsonToSsc(JsonParser.parseString(sscInJSON.substring(
                          sscInJSON.toString()
                                   .indexOf("{")))
                                                                    .getAsJsonObject()
                                                                    .toString()));
              }
          });
        Converter.convertDocumentsToSSCs(callables, convertedSSCs, nThreads);
        System.out.println("all converted");

        return new ArrayList<>(convertedSSCs);
    }

    /**
     * Converts documents to SSCs in parallel and adds them to the given collection.
     *
     * @param callables
     * @param collection
     * @param nThreads
     *
     * @throws InterruptedException
     * @see ParallelTasks#processTasks(Collection, Consumer, int)
     * @see Collection#add(Object)
     */
    public static void convertDocumentsToSSCs(final Collection<Callable<SSC>> callables,
                                              final Collection<SSC> collection,
                                              final int nThreads) throws InterruptedException {
        ParallelTasks.processTasks(callables, collection::add, nThreads);
    }

    /**
     * Converts SSCs to documents in parallel and adds them to the given collection.
     *
     * @param callables
     * @param collection
     * @param nThreads
     *
     * @throws InterruptedException
     * @see ParallelTasks#processTasks(Collection, Consumer, int)
     * @see Collection#add(Object)
     */
    public static void convertSSCsToDocuments(final List<Callable<Document>> callables, final List<Document> collection,
                                              final int nThreads) throws InterruptedException {
        ParallelTasks.processTasks(callables, collection::add, nThreads);
    }

    public static boolean sscListToJSONFile(final List<SSC> sscList, final String pathToJsonFile) {
        try {
            System.out.println("started converting... sscListToJSONFile");
            final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJsonFile));
            bw.append("{");
            bw.newLine();
            bw.flush();

            String json;
            long sscCounter = 0;
            for (final SSC ssc : sscList) {
                json = Converter.sscToJson(ssc, sscCounter);
                bw.append(json, 1, json.length()
                        - 1);
                if (sscCounter
                        < sscList.size()
                        - 1) {
                    bw.append(",");
                }
                bw.newLine();
                bw.flush();

                sscCounter++;
            }

            bw.append("}");
            bw.flush();
            bw.close();

            System.out.println("all converted");

            return true;
        } catch (final IOException e) {
            e.printStackTrace();
        }

        return false;
    }

    public static List<Query> parseQueriesFile(final String pathToQueriesFile) {
        final List<Query> queries = new ArrayList<>();
        try {
            final BufferedReader br = new BufferedReader(new FileReader(pathToQueriesFile));
            String molecularFormula = null;
            String solvent;
            Spectrum querySpectrum = new Spectrum();
            final Iterator<String> it = br.lines()
                                          .iterator();
            String line;
            String[] signalProperties;
            int metaInfoIndex;
            // process every query spectrum in queries file
            while (it.hasNext()) {
                line = it.next();
                if (line.trim()
                        .startsWith("//")) {
                    // create new query spectrum with description and solvent
                    querySpectrum = new Spectrum();
                    querySpectrum.setSignals(new ArrayList<>());
                    querySpectrum.setDescription(line.split("//")[1].trim());
                    metaInfoIndex = 1;
                    while (it.hasNext()) {
                        line = it.next();
                        if (line.trim()
                                .isEmpty()) {
                            break;
                        }
                        if (line.trim()
                                .startsWith("//")) {
                            if (metaInfoIndex
                                    == 1) {
                                try {
                                    solvent = line.split("//")[1].trim();
                                } catch (final Exception e) {
                                    solvent = null;
                                }
                                querySpectrum.setSolvent(solvent);
                            } else if (metaInfoIndex
                                    == 2) {
                                try {
                                    molecularFormula = line.split("//")[1].trim();
                                } catch (final Exception e) {
                                    molecularFormula = null;
                                }
                            }
                            metaInfoIndex++;
                        } else {
                            signalProperties = line.trim()
                                                   .split(",");
                            if (querySpectrum.getNuclei()
                                    == null) {
                                querySpectrum.setNuclei(new String[]{signalProperties[0].trim()});
                            }
                            querySpectrum.addSignal(new Signal(querySpectrum.getNuclei(), new Double[]{
                                    Double.parseDouble(signalProperties[1].trim())}, signalProperties[2].trim(),
                                                               "signal", Double.parseDouble(signalProperties[3].trim()),
                                                               Integer.parseInt(signalProperties[4].trim()), 0));
                        }
                    }
                }

                queries.add(new Query(querySpectrum, molecularFormula));
            }
            br.close();
        } catch (final IOException e) {
            e.printStackTrace();
        }

        return queries;
    }

    public static boolean hoseCodeShiftStatisticsToJSONFile(final Map<String, Double[]> hoseCodeShifts,
                                                            final String pathToJsonFile) {
        try {
            System.out.println("started converting... hoseCodeShiftsToJSONFile");
            final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJsonFile));
            bw.append("{");
            bw.newLine();
            bw.flush();

            Document subDocument;
            String json;
            long counter = 0;
            for (final Map.Entry<String, Double[]> entry : hoseCodeShifts.entrySet()) {
                subDocument = new Document();
                subDocument.append("HOSECode", entry.getKey());
                subDocument.append("values", GSON.toJson(entry.getValue()));
                json = Converter.documentToJson(new Document(String.valueOf(counter), subDocument));
                bw.append(json, 1, json.length()
                        - 1);
                if (counter
                        < hoseCodeShifts.size()
                        - 1) {
                    bw.append(",");
                }
                bw.newLine();
                bw.flush();

                counter++;
            }

            bw.append("}");
            bw.flush();
            bw.close();

            System.out.println("all converted");

            return true;
        } catch (final IOException e) {
            e.printStackTrace();
        }

        return false;
    }

    public static Map<String, Double[]> jsonFileToHOSECodeShiftStatistics(
            final String pathToJsonFile) throws FileNotFoundException {
        return Converter.jsonToHOSECodeShiftStatistics(new BufferedReader(new FileReader(pathToJsonFile)));
    }

    private static Map<String, Double[]> jsonToHOSECodeShiftStatistics(final BufferedReader br) {
        System.out.println("started converting... jsonToHOSECodeShiftStatistics");
        final Map<String, Double[]> hoseCodeShiftStatistics = new HashMap<>();
        // add all task to do
        br.lines()
          .forEach(line -> {
              if ((line.trim()
                       .length()
                      > 1)
                      || (!line.trim()
                               .startsWith("{")
                      && !line.trim()
                              .endsWith("}"))) {
                  final StringBuilder hoseCodeShiftsStatisticInJSON = new StringBuilder();
                  if (line.endsWith(",")) {
                      hoseCodeShiftsStatisticInJSON.append(line, 0, line.length()
                              - 1);
                  } else {
                      hoseCodeShiftsStatisticInJSON.append(line);
                  }
                  final JsonObject jsonObject = JsonParser.parseString(hoseCodeShiftsStatisticInJSON.substring(
                          hoseCodeShiftsStatisticInJSON.toString()
                                                       .indexOf("{")))
                                                          .getAsJsonObject();
                  hoseCodeShiftStatistics.put(jsonObject.get("HOSECode")
                                                        .getAsString(), GSON.fromJson(jsonObject.get("values")
                                                                                                .getAsString(),
                                                                                      new TypeToken<Double[]>() {
                                                                                      }.getType()));
              }
          });
        System.out.println("all converted");

        return hoseCodeShiftStatistics;
    }
}
