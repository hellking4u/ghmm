#module that reads HMMER file and converts it to XML 

import sys,re,string,StringIO
from xml.dom import minidom

def gotoLine(f,res):
    rematch = None
    r = " "
    while (rematch is None) and (r != ""):
        r = f.readline()
        rematch = res.search(r)
    if ( rematch == None ):
        print "Fehler in gotoLine"
        # sys.exit(1)
    else:
        return rematch

def build_matrix(n,m):
    "builds n x m matrix with lists"
    matrix = range(n)
    for i in range(n):
        matrix[i] = range(m)
        for j in range(m):
            matrix[i][j] = 0
    return matrix

def sum_mrows(mat):
    "sums the rows of a matrix"
    mout = range(len(mat))
    for i in range(len(mat)):
        s = 0
        for j in range(len(mat[i])):
            s=s+mat[i][j]
        mout[i] = s
    return mout

def norm_mat(mat):
    #print mat
    for i in range(len(mat)):
        s = 0.0
        #print "i: " + str(i)
        for j in range(len(mat[i])):
            s = s+mat[i][j]
        for j in range(len(mat[i])):
            #print "j: " + str(j) + ",  "+ str(mat[i][j])
            mat[i][j] = mat[i][j]/s

def red_mat_end(mat,r):
    "delete <r> rows and columns from the end of the matrix"
    for i in range(len(mat)-r):
        del mat[i][-1*r:]
    del mat[-1*r:]

def del_mat(mat,r):
    "deletes the <r>th columns and row from the matrix"
    for i in range(len(mat)):
        del mat[i][r]
    del mat[r]

def toint(i):
    "return integer if i is value or 0"
    try:
        iout = int(i)
    except ValueError:
        iout = 0
    return iout

def map_entries(dic,lis):
    "translates the letters to the number of the columns"
    dicout = {}
    for k in dic.keys():
        dicout[k] = []
        for i in range(len(dic[k])):
            dicout[k].append((lis.index(dic[k][i][0]),lis.index(dic[k][i][1])))
    return dicout

def xml_newdatanode(doc,nodename,attributename,attribute,text):
    "returns node with attribute and text"
    nod = doc.createElement(nodename)
    nod.setAttribute(attributename,attribute)
    nod.appendChild(doc.createTextNode(text))
    return nod

def write_file(strf,strcontent):
    "writes <strcontent> in file <strf>"
    try:
        f = open(strf,"w")
    except IOError,info:
        sys.stderr.write(str(info) + "\n")
        sys.exit(1)
    try:
        f.write(strcontent)
    finally:
        f.close()

class hmmer:
    "reads hmmer file and converts it to xml"

    #a few constants
    intscale = 1000.0 #hmmer format specific (but it has to be a float!!!)
    
    #positions of the entries
    lisHead = ["N","B","E","J","C","T","M","I","D"]
    lisMID = ["M","I","D"]

    #map table for HMMER file format
    dicTEntries = { "XT": [("N","B"),("N","N"),("E","C"),("E","J"),
                           ("C","T"),("C","C"),("J","B"),("J","J")],
                    "BD": [("B","D")],
                    "HMM":[("M","M"),("M","I"),("M","D"),("I","M"),
                            ("I","I"),("D","M"),("D","D"),("B","M"),("M","E")],
                    }

    #dictionary for the letters
    dicLetters = { 4:  ["A","C","G","T"],
                   20: ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q",
                        "R","S","T","V","W","Y"]}
    
    #map table translated with positions
    dicTE = map_entries(dicTEntries,lisHead)

    dicType = {"Amino":20,"Nucleotide":4}

    #coordinate offset for graphical presentation
    off_distx = 75
    off_disty = 75

    def __init__(self,strf):
        #print "start __init__"
        try:
		    f = StringIO.StringIO(strf)
        except IOError,info:
                 sys.stderr.write(str(info) + "\n")
            #     sys.exit(1)
        try:
            #print "00000"
            r = f.readline()
			#get number of match states
            self.acc  = gotoLine(f,re.compile("^ACC\s+(\w+)" ) ).group(1)
            n = int(gotoLine(f,re.compile("^LENG\s*(\d+)")).group(1))
            self.n = n
            #get type of profile hmm: amino/nucleotide
            m = self.dicType[gotoLine(f,re.compile("^ALPH\s*(\S+)")).group(1)];
            self.m = m
            #build matrix for transition: N B E J C T M1 I1 D1 M2 I2 D2 ... Mn In Dn
            self.matrans = build_matrix(3*n+6,3*n+6)
            
			#emission matrix: match state, insert state, null_model
            self.maems = [build_matrix(n,m),build_matrix(n,m),build_matrix(1,m)]
            #print "11111"  			
			#get line "XT" transitions
            trans = string.split(gotoLine(f,re.compile("^XT([\s\d\S]*)")).group(1))
            self.set_matrix("XT",trans)
            #null model
            trans = string.split(gotoLine(f,re.compile("^NULT([\s\d\S]*)")).group(1))
            self.manull = [self.H2P(trans[0],1.0),self.H2P(trans[1],1.0)] #[G->G,G->F]
            trans = string.split(gotoLine(f,re.compile("^NULE([\s\d\S]*)")).group(1))
            #print "22222"
            for k in range(len(trans)):
                self.maems[2][0][k] = self.H2P(trans[k],1/float(m))
            #main section
            gotoLine(f,re.compile("^HMM"))
            f.readline()
            #print "33333"
            #get B->D transition probability
            self.set_matrix("BD",[string.split(f.readline())[2]])
            #recurse over number of match states
            for i in range(n):
                for j in range(3):
                    lis = string.split(f.readline())
                    del lis[0]
                    if j==2:
                        #state transition probs
                        self.set_matrix("HMM",lis,i*3)
                    else:
                        #emmission probs: [match=0/insert=1][state][emission-letter]
                        for k in range(self.m):
                            #print "state: "+ str(i) +" symbol: " +str(k) +" - score = " + str(lis[k]) + ", NP = " + str(self.maems[2][0][k]) +" -> " + str(self.H2P(lis[k],self.maems[2][0][k]))							
                            self.maems[j][i][k] = self.H2P(lis[k],self.maems[2][0][k])
                            #print "null prob: " + str(self.maems[2][0][k])
            #print "44444"
            del self.maems[1][-1] #delete last (not existing) insertion state
            self.matrans[self.lisHead.index("T")][self.lisHead.index("T")] = 1.0 #set end state as loop to prevent normalization issues
            del_mat(self.matrans,-2)
            self.matrans[-1][self.lisHead.index("E")] = 1.0 # set last D->E=1
        finally:
            #pass
            f.close()
		
        #print "55555"	
        #normalize matrices:
        for i in range(len(self.maems)):
            #print "index1: " +str(i)
            #print self.maems[i]
            norm_mat(self.maems[i])
        #print "trans:"
        norm_mat(self.matrans)
		
        #print "66666"
        
        self.matrans[self.lisHead.index("T")][self.lisHead.index("T")] = 0.0
        #print "parsen fertig"
		
    def __str__(self):
        print "oben"
        hmm_str = "N= " + str(self.n)  +", M= " + str(self.m) + "\n"
        hmm_str += "Transitions: \n"
        for row in self.matrans:
            hmm_str += str(row) + "\n" 
        hmm_str += "Emissions: \n"		
        for row in self.maems:
            hmm_str += str(row) + "\n" 
        print hmm_str	
        return hmm_str


    def set_matrix(self,type,lis,offset=0):
        "fills matrix with values from line <type> and adds optional offset"
        for k in range(len(lis)):
            if lis[k]!="*":
                x_c = self.dicTEntries[type][k][0]
                y_c = self.dicTEntries[type][k][1]
                x = self.dicTE[type][k][0]+offset*(x_c in self.lisMID)
                y = self.dicTE[type][k][1]
                if y_c in self.lisMID:
                    y = y+offset
                    if ((x_c in ["D","I"]) or (x_c==y_c=="M") or (x_c=="M" and y_c=="D")) and not (x_c==y_c=="I"):
                        y=y+3
                self.matrans[x][y] = self.H2P(lis[k],1.0)
               # print (x_c + str(offset),y_c + str(offset)),self.matrans[x][y]

    def H2P(self,score,null_prob):
        "returns the probability"
        if score == "*":
            return 0
        else:
            return null_prob * 2**(float(score)/self.intscale)

    def get_dom(self):
        "returns DOM object"
        doc = minidom.Document()
        nodgraphml = doc.createElement("graphml")

        nodkey = doc.createElement("key")
        nodkey.setAttribute("for","node")
        nodkey.setAttribute("gd:type","HigherDiscreteProbDist")
        nodkey.setAttribute("id","emissions")
        nodgraphml.appendChild(nodkey)

        #declare dummy class
        nodclass = doc.createElement("hmm:class")
        nodclass.setAttribute("hmm:high","0")
        nodclass.setAttribute("hmm:low","0")
        nodmap = doc.createElement("map")
        nodsym = doc.createElement("symbol")
        nodsym.setAttribute("code","0")
        nodsym.setAttribute("desc","Simple")
        nodsym.appendChild(doc.createTextNode("N"))
        nodmap.appendChild(nodsym)
        nodclass.appendChild(nodmap)
        nodgraphml.appendChild(nodclass)

        #declare alphabet
        nodalpha = doc.createElement("hmm:alphabet")
        nodalpha.setAttribute("hmm:high",str(self.m-1))
        nodalpha.setAttribute("hmm:low","0")
        nodalpha.setAttribute("hmm:type","discrete")
        nodmap = doc.createElement("map")
        for k in self.dicLetters[self.m]:
            nodmap.appendChild(xml_newdatanode(doc,"symbol","code",str(self.dicLetters[self.m].index(k)),str(k)))
        nodalpha.appendChild(nodmap)
        nodgraphml.appendChild(nodalpha)

        dicHash = {}
        #nodes/states
        nodgraph = doc.createElement("graph")
        #add match/insert/delete states
        for k in self.lisHead:
            if k in self.lisMID:
                n = self.n
                if k == "I":
                    n=n-1
            else:
                n = 1
            for i in range(n):
                pos = self.lisHead.index(k)
                nodnode = doc.createElement("node")
                #node:data:label
                strlabel = k
                if (n>1) or (k in self.lisMID):
                    strlabel = strlabel + str(i)
                nodnode.setAttribute("id",strlabel)
                h = pos + (i*3)
                if (k=="D") and (i==n-1):
                    h=h-1 #because last I column is deleted and therefore delete column index is old one minus one
                dicHash[h] = strlabel
                nodnode.appendChild(xml_newdatanode(doc,"data","key","label",strlabel))
                #node:data:class
                nodnode.appendChild(xml_newdatanode(doc,"data","key","class","0"))
                #node:data:initial
                nodnode.appendChild(xml_newdatanode(doc,"data","key","initial",str(k=="N")))
                #node:data:ngeom
                noddata = doc.createElement("data")
                noddata.setAttribute("key","ngeom")
                nodpos = doc.createElement("pos")
                if k == "N":
                    strx = self.off_distx
                    stry = 4*self.off_disty
                elif k=="B":
                    strx = 2*self.off_distx
                    stry = 4*self.off_disty
                elif k=="J":
                    strx = (self.n+4)/2*self.off_distx
                    stry = 7*self.off_disty
                elif k=="E":
                    strx = (self.n + 3)*self.off_distx
                    stry = 2*self.off_disty
                elif k=="C":
                    strx = (self.n + 4)*self.off_distx
                    stry = 2*self.off_disty
                elif k=="T":
                    strx = (self.n + 5)*self.off_distx
                    stry = 2*self.off_disty
                elif k=="M":
                    strx = (3+i)*self.off_distx
                    stry = 3*self.off_disty
                elif k=="I":
                    strx = (3+i)*self.off_distx
                    stry = 1*self.off_disty
                elif k=="D":
                    strx = (3+i)*self.off_distx
                    stry = 5*self.off_disty
                nodpos.setAttribute("x",str(strx))
                nodpos.setAttribute("y",str(stry))
                noddata.appendChild(nodpos)
                nodnode.appendChild(noddata)
                #node:data:emissions
                if k in ["M","I"]:
                    strem = string.join(map(str,self.maems[self.lisMID.index(k)][i]),", ")
                elif k in ["N","C"]:
                    strem = string.join(map(str,self.maems[2][0]),", ")
                else:
                    strem = (self.m-1)*"0.0, " + "0.0"
                nodnode.appendChild(xml_newdatanode(doc,"data","key","emissions",strem))
                nodgraph.appendChild(nodnode)
                
        #edges/transitions
        for i in range(len(self.matrans)):
            for j in range(len(self.matrans[i])):
                if self.matrans[i][j]!=0:
                    nodedge = doc.createElement("edge")
                    nodedge.setAttribute("source",dicHash[i])
                    nodedge.setAttribute("target",dicHash[j])
                    nodedge.appendChild(xml_newdatanode(doc,"data","key","prob",str(self.matrans[i][j])))
                    nodgraph.appendChild(nodedge)

        nodgraphml.appendChild(nodgraph)

        doc.appendChild(nodgraphml)
        return doc

    def write_prettyxml(self,strf):
        "writes hmm in pretty xml format to file"
        write_file(strf,self.get_dom().toprettyxml())

    def write_xml(self,strf):
        "writes hmm in xml format to file"
        write_file(strf,self.get_dom().toxml())
