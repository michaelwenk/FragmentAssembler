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
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
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
        document.append("rootAtomIndex", ssc.getRootAtomIndex());
        document.append("hoseCodes", ssc.getHoseCodes());
        document.append("unsaturatedAtomIndices", ssc.getUnsaturatedAtomIndices());

        return document;
    }


    public static SSC jsonToSsc(final String json) {
        final JsonObject jsonObject = new JsonParser().parse(json)
                                                      .getAsJsonObject();
        return new SSC(GSON.fromJson(jsonObject.get("structure"), ExtendedConnectionMatrix.class)
                           .toAtomContainer(), GSON.fromJson(jsonObject.get("spectrum"), Spectrum.class),
                       GSON.fromJson(jsonObject.get("assignment"), Assignment.class), jsonObject.get("rootAtomIndex")
                                                                                                .getAsInt(),
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
        System.out.println("started converting...");

        final ConcurrentLinkedQueue<SSC> convertedSSCs = new ConcurrentLinkedQueue<>();
        final JsonParser jsonParser = new JsonParser();
        //        final JsonReader jsonReader = GSON.newJsonReader(br);
        //        jsonReader.setLenient(true);
        //        System.out.println(jsonParser.parse(jsonReader)
        //                                     .getAsJsonObject());
        final List<Callable<SSC>> callables = new ArrayList<>();
        // add all task to do
        br.lines()
          .parallel()
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
                  callables.add(() -> Converter.jsonToSsc(jsonParser.parse(sscInJSON.substring(sscInJSON.toString()
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

    /**
     * Stores this SSC library into a string in JSON format.
     *
     * @return this SSC library as JSON string
     */
    public static String sscListToJson(final List<SSC> sscList) {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("{");
        stringBuilder.append("\n");
        String json;
        long sscCounter = 0;
        for (final SSC ssc : sscList) {
            json = Converter.sscToJson(ssc, sscCounter);
            stringBuilder.append(json, 1, json.length()
                    - 1);
            if (sscCounter
                    < sscList.size()
                    - 1) {
                stringBuilder.append(",");
            }
            stringBuilder.append("\n");

            sscCounter++;
        }
        stringBuilder.append("}");

        return stringBuilder.toString();
    }

    public static boolean sscListToJSONFile(final List<SSC> sscList, final String pathToJsonFile) {
        try {
            final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJsonFile));
            bw.write(Converter.sscListToJson(sscList));
            bw.close();

            return true;
        } catch (final IOException e) {
            e.printStackTrace();
            return false;
        }
    }

    public static List<Query> parseQueriesFile(final String pathToQueriesFile) {
        final List<Query> queries = new ArrayList<>();
        try {
            final BufferedReader br = new BufferedReader(new FileReader(pathToQueriesFile));
            String molecularFormula = "";
            Spectrum querySpectrum = new Spectrum();
            final Iterator<String> it = br.lines()
                                          .iterator();
            String line;
            String[] signalProperties;
            // process every query spectrum in queries file
            while (it.hasNext()) {
                line = it.next();
                if (line.trim()
                        .startsWith("//")) {
                    // create new query spectrum with description
                    querySpectrum = new Spectrum();
                    querySpectrum.setSignals(new ArrayList<>());
                    querySpectrum.setDescription(line.split("//")[1].trim());
                    molecularFormula = "";
                    while (it.hasNext()) {
                        line = it.next();
                        if (line.trim()
                                .isEmpty()) {
                            break;
                        }
                        if (line.trim()
                                .startsWith("//")) {
                            try {
                                molecularFormula = line.split("//")[1].trim();
                                System.out.println(" read molecular formula: "
                                                           + molecularFormula);
                            } catch (final Exception e) {
                                molecularFormula = null;
                            }
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
}
