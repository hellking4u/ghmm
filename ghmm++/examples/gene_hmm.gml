<graphml>

 <desc>Class HMM Genefinder</desc>
 <key id="label" for="node"><desc>A human readable description of the vertex</desc></key>
 <key id="initial" hmm:type="float" for="node"><desc>Initial probability of a state</desc></key>
 <key id="order" hmm:type="int" for="node"><desc>Order of a state</desc></key>
 <key id="class" hmm:type="hmm:class" for="node"><desc>Class of a state</desc></key>
 <key id="tiedto" hmm:type="hmm:id" for="node"><desc>Which state is it tied to</desc></key>
 <key id="emissions" gd:type="HigherDiscreteProbDist" for="node"><desc>Discrete emission probability</desc></key>
 <key id="prob" hmm:type="float" for="edge"><desc>transition probability</desc></key>

 <!-- ====== Defaults for the visual representation of nodes/edges -->

 <!-- node attribute of type paint with default value -->
 <key id="npaint" for="node" gd:type="paint">
  <paint red="255" blue="228" green="255"/>
 </key>

 <!-- edge attribute of type paint with default value -->
 <key id="epaint" for="edge" gd:type="paint">
  <paint red="0" blue="0" green="0" style="solid"/>
 </key>

 <!-- node attribute of type point with default value -->
 <key id="ngeom" for="node" gd:type="point">
  <point shape="circle" width="25" height="25"/>
 </key>

 <!-- edge attribute of type line with default value -->
 <key id="egeom" for="edge" gd:type="line">
  <line shape="poly" width="1"/>
 </key>

 <hmm:class hmm:low="0" hmm:high="4"> 
   <map>
     <symbol code="0" desc="Intergenic">N</symbol>
     <symbol code="1" desc="Single">S</symbol>
     <symbol code="2" desc="First">F</symbol>
     <symbol code="3" desc="Internal">I</symbol>
     <symbol code="4" desc="Terminal">T</symbol>
   </map>
 </hmm:class>

 <hmm:alphabet hmm:type="discrete" hmm:low="0" hmm:high="4"> 
   <map>
     <symbol code="0">A</symbol>
     <symbol code="1">C</symbol>
     <symbol code="2">G</symbol>
     <symbol code="3">T</symbol>
     <symbol code="4">N</symbol>
   </map>
 </hmm:alphabet>

 <!-- ====== The HMM ========================================== -->
 <graph>
     <node id="1">
       <data key="label">Intron</data>
       <data key="class">0</data>
       <data key="initial">0.7</data>
       <data key="order">0</data>
       <data key="emissions">0.2, 0.2, 0.5, 0.1, 0.0</data>
       <data key="ngeom"><pos x="100" y="100"/></data>
     </node>
     <node id="2">
       <data key="label">Exon</data>
       <data key="class">1</data>
       <data key="initial">0.3</data>
       <data key="order">0</data>
       <data key="ngeom"><pos x="200" y="200"/></data>
       <data key="emissions">0.9, 0.1, 0.0, 0.0, 0.0</data>
     </node>
     <edge source="1" target="2">
       <data key="prob">0.5</data>
     </edge>
     <edge source="1" target="1">
       <data key="prob">0.5</data>
     </edge>
     <edge source="2" target="2">
       <data key="prob">0.2</data>
     </edge>
     <edge source="2" target="1">
       <data key="prob">0.8</data>
     </edge>
</graph>
</graphml>

