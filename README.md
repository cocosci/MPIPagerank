# MPJ-PagRank

javac -cp .:$MPJ_HOME/lib/mpj.jar MPIPageRank.java 

mpjrun.sh -np 4 MPIPageRank input.txt out.txt 1000 0.85



# Result

```
[zhenk@silo src]$ mpjrun.sh -np 4 MPIPageRank input.txt out.txt 1000 0.85 
Picked up _JAVA_OPTIONS: -Xms512m -Xmx512m
MPJ Express (0.44) is started in the multicore configuration
Picked up _JAVA_OPTIONS: -Xms512m -Xmx512m
For Process Rank 0:
Running time for loading for :13ms
Running time for calculation:1524ms
For Process Rank 2:
Running time for loading for :119msz`
Running time for calculation:1418ms
For Process Rank 3:
Running time for loading for :118ms
Running time for calculation:1418ms
For Process Rank 1:
Running time for loading for :116ms
Running time for calculation:1417ms
1 :-  0.38440094881355413
2 :-  0.34291028550837954
4 :-  0.08088569323449774
3 :-  0.039087092099966095
5 :-  0.039087092099966095
0 :-  0.03278149315934399
6 :-  0.016169479016858404
7 :-  0.016169479016858404
8 :-  0.016169479016858404
9 :-  0.016169479016858404
```
