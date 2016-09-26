import mpi.MPI;

import java.io.*;
import java.util.*;
import java.util.stream.*;
import java.nio.file.*;
import java.nio.charset.Charset;

public class MPIPageRank {
    // adjacency matrix read from file
    private HashMap<Integer, ArrayList<Integer>> adjMatrix = new HashMap<Integer, ArrayList<Integer>>();
    // input file name
    private String inputFile = "";
    // output file name
    private String outputFile = "";
    // number of iterations
    private int iterations = 10;
    // damping factor
    private double df = 0.85;
    // number of precess, rank and size
    private int numberofProcess = 1;
    private int rank = 0;
    private int size = 1;
    private int numofLine = 0;

    private double[] localPageRankValues,globalPageRankValue, danglingPageRankValue;



    /**
     * Parse the command line arguments and update the instance variables. Command line arguments are of the form
     * <input_file_name> <output_file_name> <num_iters> <damp_factor>
     *
     * @param args arguments
     */
    public void parseArgs(String[] args) {
        inputFile=args[3];
        outputFile=args[4];
        iterations=Integer.parseInt(args[5]);
        df=Double.parseDouble(args[6]);
        numberofProcess = Integer.parseInt(args[1]);
        System.out.print(numberofProcess);
        // the fifth parameter deals with the number of processes.
        MPI.Init(args);
        rank = MPI.COMM_WORLD.Rank();
        size = MPI.COMM_WORLD.Size();


    }

    public int getrank() {
        return rank;
    }

    public void setFileLength() throws IOException {
        try (Stream<String> lines = Files.lines(Paths.get(inputFile), Charset.defaultCharset())) {
            numofLine = (int)lines.count();
        }

    }

    public void initializePageRankVariables(){
        globalPageRankValue=new double[numofLine];
        danglingPageRankValue=new double[numofLine];
        localPageRankValues=new double[numofLine];
        for (int i = 0; i < numofLine; i++) {
            double thesize = numofLine;
            globalPageRankValue[i]= 1 / thesize;
            danglingPageRankValue[i] = 0;
            localPageRankValues[i] = 0;
        }
    }

    public int getChunkSize(){
        return numofLine % numberofProcess == 0 ? numofLine / numberofProcess : numofLine / numberofProcess + 1;
    }

    public void inititalizeAdjMatrix(ArrayList<String> receivedChunk){
        for (int i2 = 0; i2 < receivedChunk.size(); i2++) {
            String line = receivedChunk.get(i2);
            String split_line[] = line.split(" ");
            if (split_line.length > 1) {
                ArrayList<Integer> nodes = new ArrayList<>();
                for (int i = 1; i < split_line.length; i++)
                    nodes.add(Integer.parseInt(split_line[i]));
                adjMatrix.put(Integer.parseInt(split_line[0]), nodes);
            } else if (split_line.length == 1)
                adjMatrix.put(Integer.parseInt(split_line[0]), null);
        }
    }
    /**
     * Read the input from the file and populate the adjacency matrix
     * <p>
     * The input is of type
     * <p>
     * The first value in each line is a URL. Each value after the first value is the URLs referred by the first URL.
     * For example the page represented by the 0 URL doesn't refer any other URL. Page
     * represented by 1 refer the URL 2.
     *
     * @throws java.io.IOException if an error occurs
     */

    public void loadInput() throws IOException {

        setFileLength();
        // init two global structures.

        initializePageRankVariables();

        try (BufferedReader br = new BufferedReader(new FileReader(inputFile))) {
            // the remainder may not be zero!!!!
            // either read all in a global adjmat and scatter/read certain lines from file and then send the big string.
            // if 10 lines and 4 processes, then it's 3 3 3 1.
            int chunkSize = getChunkSize();
            int count = 0;
            int senderIndex = 1; //start from process rank1.
            boolean fulfilledRank0 = false;
            ArrayList<String> chunkstrrecv = new ArrayList<String>();
            Object[] chunkstrrecvobj = new Object[1];
            if (rank == 0) {
                System.out.println("Chunk Size:"+ chunkSize);

                ArrayList<String> chunkstr = new ArrayList<String>();
                for (int i = 0; i < numofLine; i++) {
                    //each time add one line from the input file to the chunkstr.
                    String tempstr = br.readLine();
                    chunkstr.add(tempstr);
                    count++;
                    //if we fill up the chunk, we send the chunksize to the required process
                    if (count == chunkSize || i == numofLine - 1) {
                        //we first init the process 0;
                        if (!fulfilledRank0) {
                            chunkstrrecv.addAll(chunkstr);
                            fulfilledRank0 = true;
                        } else {
                            //after initializing process 0, send the object of array<string> to the process with the number of sendtowho.
                            Object[] chunkstrobj = new Object[1];
                            chunkstrobj[0] = (Object) chunkstr;
                            MPI.COMM_WORLD.Send(chunkstrobj, 0, 1, MPI.OBJECT, senderIndex, 1);
                            senderIndex++;
                        }
                        chunkstr.clear();
                        count = 0;
                    }
                }
            }else
                MPI.COMM_WORLD.Recv(chunkstrrecvobj, 0, 1, MPI.OBJECT, 0, 1);
            if (rank != 0)
                chunkstrrecv = (ArrayList<String>) chunkstrrecvobj[0];
            inititalizeAdjMatrix(chunkstrrecv);
        }
    }

    /**
     * Do fixed number of iterations and calculate the page rank values. You may keep the
     * intermediate page rank values in a hash table.
     */

    // for each process, it's roughly the same.
    public void calculatePageRank() {

        for (int i = 0; i < iterations; i++) {
            for (int j = 0; j < numofLine; j++)
                localPageRankValues[j] = 0.0;
            double danglingContrib = 0.0;
            Iterator it = adjMatrix.entrySet().iterator();
            while (it.hasNext()) {
                Map.Entry<Integer, List> pair = (Map.Entry) it.next();

                //If it is a dangling node,
                if (pair.getValue() == null)
                    danglingContrib += globalPageRankValue[pair.getKey()]/numofLine;
                else{
                    int current_size = pair.getValue().size();
                    Iterator iter = pair.getValue().iterator();
                    // For each outbound link for a node
                    while (iter.hasNext()) {
                        int node = Integer.parseInt(iter.next().toString());
                        double temp = localPageRankValues[node];
                        temp += globalPageRankValue[pair.getKey()] / current_size;
                        localPageRankValues[node] = temp;
            }}}
            double tempSend[] = new double[1];
            double tempRecv[] = new double[1];
            tempSend[0] = danglingContrib;
            //Send the dangling Contribution
            MPI.COMM_WORLD.Allreduce(tempSend,0,tempRecv,0,1,MPI.DOUBLE,MPI.SUM);
            //Send the localPageRankValue
            MPI.COMM_WORLD.Allreduce(localPageRankValues,0,globalPageRankValue,0,numofLine,MPI.DOUBLE,MPI.SUM);
            if(rank==0){
                for(int k=0;k<numofLine;k++) {
                    globalPageRankValue[k] += tempRecv[0];
                    globalPageRankValue[k] = df * globalPageRankValue[k] + (1 - df) * (1.0 / (double) numofLine);
                }
            }
            MPI.COMM_WORLD.Bcast(globalPageRankValue, 0, numofLine, MPI.DOUBLE, 0);
        }
    }
    /**
     * Print the pagerank values. Before printing you should sort them according to decreasing order.
     * Print all the values to the output file. Print only the first 10 values to console.
     *
     * @throws IOException if an error occurs
     */

    public void printValues() throws IOException {
        for(int k=0;k<numofLine;k++)
            System.out.println(k+":"+globalPageRankValue[k]);
        HashMap<Integer,Double> pageRank= new HashMap<Integer,Double>();
        for (int i = 0; i < globalPageRankValue.length; i++) {
            pageRank.put(i, globalPageRankValue[i]);
        }

        Set<Map.Entry<Integer, Double>> set = pageRank.entrySet();


        List<Map.Entry<Integer, Double>> list = new ArrayList<Map.Entry<Integer, Double>>(set);
        Collections.sort(list, new Comparator<Map.Entry<Integer, Double>>() {
            public int compare(Map.Entry<Integer, Double> o1, Map.Entry<Integer, Double> o2) {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });
        try (Writer writer = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(outputFile), "utf-8"))) {
            int i = 1;
            Iterator it = list.iterator();
            while (it.hasNext()) {
                Map.Entry<Integer, Double> pair = (Map.Entry) it.next();
                System.out.println(pair.getKey() + " :-  " + pair.getValue());
                writer.write(pair.getKey() + " :-  " + pair.getValue() + "\n");
                i++;
                if (i == 11)
                    break;
            }
        }
    }

    /**
     * Print out the private variables of adjmat to debug.
     */
    private void printadjMat() {
        // TODO Auto-generated method stub
        // for AdjMatrix
        System.out.println("For AdjMatrix..."+rank);
        Iterator it = adjMatrix.values().iterator();
        while (it.hasNext()) {
            ArrayList<Integer> arr = (ArrayList<Integer>) it.next();
            if (arr != null) {
                for (int i = 0; i < arr.size(); i++)
                    System.out.print(arr.get(i) + " ");
                System.out.println();
            } else System.out.println("null");
        }
        MPI.Finalize();
    }

    public static void main(String[] args) throws IOException {


        MPIPageRank MPIPR = new MPIPageRank();
        if (args.length < 2)
            System.exit(1);

        MPIPR.parseArgs(args);
        MPIPR.loadInput();
        MPIPR.calculatePageRank();
        if (MPIPR.getrank() == 0) {
            MPIPR.printValues();
        }

    }


}