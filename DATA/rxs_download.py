#!/usr/bin/env python3
# Python wrapper for the Reaxys API
#
# Version: 1.04
# 
# Author:  Dr. Sebastian Radestock, Elsevier
# Author:  Dr. Alexander Riemer, Elsevier
# Date:    8 June 2016

import os, http.cookiejar, urllib, xml.dom.minidom, string, sys
import codecs
import gzip
import zipfile
import re

class _Reaxys_API(object):
    def __init__(self, logger):
        self.logger = logger

        self.url = ""
        self.headers = {'Content-type' : 'text/xml; charset="UTF-8"'}
        self.callername = ""
        self.sessionid = ""
        self.resultname = ""
        self.resultsize = ""
        self.citationset = ""
        self.citationcount = ""

    def _get_resultname(self, response_xml):
        response_dom = xml.dom.minidom.parseString(response_xml)
        # Length of response_dom.getElementsByTagName("resultsname") should always be 1.
        # Node resultsname should not conatin subnodes.
        return response_dom.getElementsByTagName("resultname")[0].childNodes[0].nodeValue

    def _get_resultsize(self, response_xml):
        response_dom = xml.dom.minidom.parseString(response_xml)
        # Length of response_dom.getElementsByTagName("resultsize") should always be 1.
        # Node resultsize should not conatin subnodes.
        return response_dom.getElementsByTagName("resultsize")[0].childNodes[0].nodeValue
    
    def _get_citationset(self, response_xml):
        response_dom = xml.dom.minidom.parseString(response_xml)
        # Length of response_dom.getElementsByTagName("citationset") should always be 1.
        # Node citationset should not conatin subnodes.          
        return response_dom.getElementsByTagName("citationset")[0].childNodes[0].nodeValue
    
    def _get_citationcount(self, response_xml):
        response_dom = xml.dom.minidom.parseString(response_xml)
        # Length of response_dom.getElementsByTagName("citationcount") should always be 1.
        # Node citationcount should not conatin subnodes.          
        return response_dom.getElementsByTagName("citationcount")[0].childNodes[0].nodeValue

    def get_facts_availability(self, response_xml, field):
        facts_availability = "0"
        response_dom = xml.dom.minidom.parseString(response_xml)
        facts = response_dom.getElementsByTagName("facts")[0]
        for fact in facts.childNodes:
            if 'name="' + field + '"' in fact.toxml():
                facts_availability = fact.childNodes[0].nodeValue.split("(")[0]
        return facts_availability

    def get_field_content(self, response_xml, field):
        field_content = []
        response_dom = xml.dom.minidom.parseString(response_xml)
        for element in response_dom.getElementsByTagName(field):
            # If node contains further sub-nodes: return full xml.
            if len(element.childNodes) > 1:
                field_content.append(element.toxml())
            # If node does not contain further sub-nodes: return node value.
            elif len(element.childNodes) == 1:
                field_content.append(element.childNodes[0].nodeValue)
        return field_content

    def connect(self, url, url_main, username, password, callername):
        self.url = url
        self.callername = callername
        cookies = http.cookiejar.CookieJar()
        self.headers['Cookie'] = ""
        
        connect_template = """<?xml version="1.0"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s">
              <statement command="connect" username="%s" password="%s"/>
            </request>
          </xf>\n"""
        payload = connect_template%(callername, username, password)

        # Header reset.
        self.headers = {'Content-type' : 'text/xml; charset="UTF-8"'}

        # ELSAPI support
        self.headers['X-ELS-APIKey'] = callername
        self.headers['Accept'] = "*/*"
        request = urllib.request.Request(self.url, payload.encode("utf-8"), self.headers)
        
        #self.logger('-----------------------\nQuery from connect:')
        #self.logger(payload)

        response = urllib.request.urlopen(request)

        response_xml = response.read()

        #self.logger('-----------------------\nResponse headers from connect:')
        #self.logger(response.info())
        #self.logger('-----------------------\nResponse from connect:')
        #self.logger(response_xml.decode("utf-8"))

        # Get sessionid.
        response_dom = xml.dom.minidom.parseString(response_xml)
        element = response_dom.getElementsByTagName("sessionid")
        self.sessionid = element[0].childNodes[0].nodeValue
        
        # Cookies are read from the response and stored in self.header
        #     which is used as a request header for subsequent requests.
        cookies.extract_cookies(response, request)

        self.headers['Cookie'] = "; ".join("%s=%s" % (cookie.name, cookie.value) for cookie in cookies)


    def disconnect(self):
        disconnect_template = """<?xml version="1.0"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s">
              <statement command="disconnect" sessionid="%s"/>
            </request>
          </xf>\n"""
        payload = disconnect_template%(self.callername, self.sessionid)

        request = urllib.request.Request(self.url, payload.encode("utf-8"), self.headers)
        response = urllib.request.urlopen(request)
        response_xml = response.read()
        
        #self.logger('-----------------------\nQuery headers from disconnect:')
        #self.logger(self.headers)
        #self.logger('-----------------------\nQuery from disconnect:')
        #self.logger(payload)
        #self.logger('-----------------------\nResponse headers from disconnect:')
        #self.logger(response.info())
        #self.logger('-----------------------\nResponse from disconnect:')
        #self.logger(response_xml.decode("utf-8"))
            
    def select(self, dbname, context, where_clause, order_by, options):
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s" sessionid="">
              <statement command="select"/>
              <select_list>
                <select_item/>
              </select_list>
              <from_clause dbname="%s" context="%s">
              </from_clause>
              <where_clause>%s</where_clause>
              <order_by_clause>%s</order_by_clause>
              <options>%s</options>
            </request>
          </xf>\n"""
        payload = select_template%(self.callername, dbname, context, where_clause, order_by, options)
        
        request = urllib.request.Request(self.url, payload.encode("utf-8"), self.headers)
        response = urllib.request.urlopen(request)
        response_xml = response.read()
        
        #self.logger('-----------------------\nQuery headers from select:')
        #self.logger(self.headers)
        #self.logger('-----------------------\nQuery from select:')
        #self.logger(payload)
        #self.logger('-----------------------\nResponse headers from select:')
        #self.logger(response.info())
        #self.logger('-----------------------\nResponse from select:')
        #self.logger(response_xml.decode("utf-8"))

        self.resultname = self._get_resultname(response_xml)
        self.resultsize = self._get_resultsize(response_xml)
        
        if ("NO_CORESULT" not in options) and ("C" not in context):
            self.citationset = self._get_citationset(response_xml)
            self.citationcount = self._get_citationcount(response_xml)


    def expand(self, dbname, first_item, last_item, where_clause):
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s" sessionid="%s">
              <statement command="expand"/>
              <from_clause dbname="%s" first_item="%s" last_item="%s">
              </from_clause>
              <where_clause>%s</where_clause>
            </request>
          </xf>\n"""
        payload = select_template%(self.callername, self.sessionid, dbname, first_item, last_item, where_clause)
        
        request = urllib.request.Request(self.url, payload.encode("utf-8"), self.headers)
        response = urllib.request.urlopen(request)
        response_xml = response.read()
        
        #self.logger('-----------------------\nQuery headers from expand:')
        #self.logger(self.headers)
        #self.logger('-----------------------\nQuery from expand:')
        #self.logger(payload)
        #self.logger('-----------------------\nResponse headers from expand:')
        #self.logger(response.info())
        ##self.logger('-----------------------\nResponse from expand:')
        ##self.logger(response_xml)

        return response_xml

    def post(self, payload):
        request = urllib.request.Request(self.url, payload.encode("utf-8"), self.headers)
        response = urllib.request.urlopen(request)
        response_xml = response.read()
        
        #self.logger('-----------------------\nQuery headers from post:')
        #self.logger(self.headers)
        #self.logger('-----------------------\nQuery from post:')
        #self.logger(payload)
        #self.logger('-----------------------\nResponse headers from post:')
        #self.logger(response.info())
        ##self.logger('-----------------------\nResponse from post:')
        ##self.logger(response_xml)



    def retrieve(self, resultname, select_items, first_item, last_item, order_by, group_by, group_item, options, checkRXD):
        # if group_by is given, please provide group_item, e.g. "1" or "1,2"
        
        if group_by != '':
            grouplist = 'grouplist="' + group_item + '"'
        else:
            grouplist = ''

        select_item_template = """                <select_item>%s</select_item>\n"""
        select_template = """<?xml version="1.0" encoding="UTF-8"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s" sessionid="%s">
              <statement command="select"/>
              <select_list>\n"""
        for index in range (0,len(select_items)):
            select_template = select_template + select_item_template%(select_items[index])
        select_template = select_template + """              </select_list>
              <from_clause resultname="%s" %s first_item="%s" last_item="%s">
              </from_clause>
              <order_by_clause>%s</order_by_clause>
              <group_by_clause>%s</group_by_clause>
              <options>%s</options>
            </request>
          </xf>\n"""
        payload = select_template%(self.callername, self.sessionid, resultname, grouplist, first_item, last_item, order_by, group_by, options)
        
        request = urllib.request.Request(self.url, payload.encode("utf-8"), self.headers)
        response = urllib.request.urlopen(request)
        response_xml = response.read()

        ###
        # Check for 100 RXDs (RX IDs with more than 100 RXDs associated)
        if checkRXD:
            #t0=time()
            for rx in str(response_xml).split("<RX>")[1:]:  # Split in reactions and chop off header
                if rx.count("<RXD>") == 100:
                    rxid=re.sub(r'.*<RX.ID.*>(.+)</RX.ID>.*', r'\1', rx[:50])  # Get reaction ID to retrieve entry iteratively later.
                    self.logger("#Reiterate{}#".format(rxid))
            #print("TIME CHECKING\t{:.5f}".format(time()-t0))  # Program spends tipically 0.01s doing this. It's worth it!
        ###
        
        #self.logger('-----------------------\nQuery headers from retrieve:')
        #self.logger(self.headers)
        #self.logger('-----------------------\nQuery from retrieve:')
        #self.logger(payload)
        #self.logger('-----------------------\nResponse headers from retrieve:')
        #self.logger(response.info())
        #self.logger('-----------------------\nResponse from retrieve:')
        #self.logger(response_xml)

        return response_xml


class Downloader(object):
    def __init__(self, username, password, callername, logName, debug):
        self.username = username
        self.password = password
        self.callername = callername
        self.logName = logName
        self.debug = debug

        # Reaxys API server
        self.url_main = "www.reaxys.com"
        # Reaxys CERT API server
        self.url = "https://" + self.url_main + "/reaxys/api"

    def __enter__(self):
        def logger(s, self=self):
            with open(self.logName, "a") as f:
                f.write(str(s)+"\n")

        self.ra = _Reaxys_API(logger)
        return self._Handler(self)  # Use in a with .. as .. statement will generate Handler object, and close it after use.
    
    def __exit__(self, type, value, traceback):
        if self.debug:  self.log.close()
        else: pass

    class _Handler(object):
        class ResultChunk(object):
            def __init__(self, idStart, numEntries, xmlBytes):
                self.idStart = idStart
                self.numEntries = numEntries
                self.xmlBytes = xmlBytes

        def __init__(self, d):
            self.d = d

        def logger(self, s):
            with open(self.d.logName, "a") as f:
                f.write(str(s)+"\n")

        def _fetch(self, lFromId, lToId, context, idField, fields, check_str=False, checkRXD=False):
            d = self.d
            selectLimit = 10000
            chunkSize = 100
            for fromId in range(lFromId, lToId, selectLimit):

                d.ra.connect(d.url, d.url_main, d.username, d.password, d.callername)
                toId = min(fromId + selectLimit, lToId)
                #print("Selection [%d, %d[" % (fromId, toId))
                # select parameters:
                # DatabaseName, context, WhereClause, OrderBy, Options
                # The DatbaseName is "RX" for Reaxys
                # (the 'between' operator is given an inclusive range)
                query = "%s between %s and %s" % (idField, fromId, toId - 1)
                if check_str: query += " AND exists(YY)"

                d.ra.select("RX", context, query,"", "WORKER,NO_CORESULT")

                for chunkStart in range(fromId, toId, chunkSize):
                    # they use 1-indexing and inclusive ranges
                    offsetFirst = chunkStart - fromId + 1 
                    offsetLast = offsetFirst + chunkSize - 1

                    response = d.ra.retrieve(d.ra.resultname, fields,
                                                     "%s" % offsetFirst, "%s" % offsetLast,
                                                     "", "", "",
                                                     "ISSUE_RXN=true, ISSUE_RCT=true, COMPRESS=true", 
                                                     checkRXD=checkRXD)

                    # yield chunk result and RX.IDs with more than 100 RXDs
                    yield self.ResultChunk(chunkStart, chunkSize, response)
                d.ra.disconnect()

        def fetchCompounds(self, fromId, toId, check_str):
            return self._fetch(fromId, toId, "S", "IDE.XRN", ["PHY","SPE","PDE","ECO","ID"], check_str)

        def fetchReactions(self, fromId, toId, n, m, checkRXD):
            return self._fetch(fromId, toId, "R", "RX.ID", ["RX","RY","RXD({},{})".format(n,m)],checkRXD=checkRXD)

        def fetchXCitations(self, fromId, toId, check_str):
            return self._fetch(fromId, toId, "C", "CNR.CNR", ["CIT"], check_str)

        def fetchBioactivity(self, fromId, toId, check_str):
            return self._fetch(fromId, toId, "DPI", "DAT.ID", ["DAT"], check_str)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Reaxys XML downloader.")

    parser.add_argument("--type","-t", choices=["S", "R", "C" , "B"], required=True,
                        help="Which type to download. S=Substances, R=Reactions, C=Citations , B=Bioactivities")
    parser.add_argument("--num","-n", type=int, required=True, 
                        help="The size of the chunck of entries to download.")
    parser.add_argument("--offset","-o", type=int, required=False, default=0,
                        help="Which id to start from.")
    parser.add_argument("--name", required=False, default="log",
                        help="The name of the archive to produce.")
    parser.add_argument("--debug", required=False, default=False,
                        help="Produce verbose output for each operation, for debugging purposes.")
    parser.add_argument("--check_str","-s", required=False, action="store_true",
                        help="Retrieve only entries that contain a substance's structure.\n \
                             Type the flag for using it, else defaults to not using the option.")
    parser.add_argument("--nproc","-N", required=False, default=1,type=int,
                        help="Number of processors to use. Up to 10 as not more logins are allowed by Reaxys.")
    parser.add_argument("--check_rxd","-x", required=False, action="store_true",
                        help="Check for reactions with possibly more than 100 RXDs.\n \
                             Type the flag for using it, else defaults to not using the option.")
    parser.add_argument("--checkpoint","-k", required=False, action="store_true",
                        help="Restart download from checkpoint. Check the files in dir and restart from there.\n \
                             Type the flag for using it, else defaults to not using the option.")
    args = parser.parse_args()

    from time import time
    t0 = time()

    def getData(proc, size, todir):

        try:
            if not os.path.isdir('{}/p{}'.format(todir,proc)):
                os.mkdir('{}/p{}'.format(todir,proc))

            logName = "{}/{}_p{}.log".format(todir, args.name, proc)
            
            if args.checkpoint:  # Check whether this process finished or not
                # Last line should read "***Finished without errors"
                with open(logName, 'r') as f:
                    for line in f: pass
                if "*** Finished without errors." in line:
                    return 0  # Return success again.

            downloader = Downloader(
                # Reaxys credentials (leave blank if IP-based authentication is used)
                username ="***REMOVED***",
                password = "***REMOVED***",
                # API Key
                callername = "uni_leipzig_dec2015",
                logName = logName,
                debug=args.debug
            )

            with downloader as d:

                # Calculate what the ID boundaries are for this processor
                ChunkSize = args.num//size
                fromId = args.offset + proc*ChunkSize
                toId = args.offset + (proc+1)*ChunkSize
                if proc==size-1:  # For last process toId is just args.end
                    toId = args.offset + args.num

                if args.checkpoint: # Check what the progress is so far in each subdir p{i}.
                    # List all elements in todir/p{i} and find maximum number so far.
                    files = os.listdir("{}/p{}/".format(todir, proc))
                    f=[i for i in files if "_c" not in i]  # Check only on non-marked files
                    if len(f)>0:
                        progress = max(list(map(lambda x: int(re.sub(r'n([0-9]+).xml', r'\1',x)), f)))
                        fromId = progress-100 # Restart from one file before just in case last one wasn't ok
                        d.logger("p{}: Restarting from ID {}".format(proc, fromId))
                    else: # If no file was written before, start from original fromId
                        d.logger("p{}: Restarting from ID {}".format(proc, fromId))


                # Proceed with retrieval
                d.logger("proc {}: fromId: {} toID {}".format(proc,fromId,toId))
                if args.type == "S":
                    res = d.fetchCompounds(fromId, toId, args.check_str)
                elif args.type == "R":
                    res = d.fetchReactions(fromId,toId,1,100,args.check_rxd)
                elif args.type == "C":
                    res = d.fetchXCitations(fromId, toId, args.check_str)
                elif args.type == "B":
                    res = d.fetchBioactivity(fromId,toId,args.check_str)    
                else:   assert False


                ####### This is main iteration: Correspond to pass through all indices
                for chunk in res:   # This is the loop calling the generator!!
                    with open("{}/p{}/n{}.xml".format(todir, proc, chunk.idStart), "wb") as f:
                        f.write(chunk.xmlBytes)

                ####### Corrections are done **after** main loop is finished! 
                ####### Here iterate over selected RX IDs 

                repeat_ids = []
                with open(logName, "r") as logf:
                    for line in logf:
                        if re.search("#Reiterate",line):
                            repId = int(re.sub(r'#Reiterate([0-9]+)#',r'\1',line))
                            repeat_ids.append(repId)

                for rxid in repeat_ids:
                    d.logger("\t***Reiterating over RX ID {}".format(rxid))                
                    iternm=True
                    n,m=1,100
                    
                    while iternm:
                        d.logger("proc {}: fromId {} toID {}".format(proc,rxid, rxid+1))
                        res = d.fetchReactions(rxid, rxid+1, n, m, checkRXD=False)

                        for chunk in res:   # This is the loop calling the generator!!
                            with open("{}/p{}/n{}_c.xml".format(todir, proc, chunk.idStart), "ab") as f:
                                f.write(chunk.xmlBytes)
                                if str(chunk.xmlBytes).count("<RXD>")<100:  # If number of <RXD> is <100, stop loop
                                    iternm=False
                                else: # Get next chunk of RXDs
                                    n+=100; m+=100
                ###############

            d.logger("*** Finished without errors.")
            return 0 # In case everything finished without any problem.

        except: return -1  # Return error status
        


    import multiprocessing as mp

    if args.type == "R":        todir="./DATA/RXN/"
    elif args.type == "S" :     todir="./DATA/SUB/"
    elif args.type == "C" :     todir="./DATA/CIT/"
    elif args.type == "B" :     todir="./DATA/DPI/"


    if not os.path.isdir('{}'.format(todir)):
        os.mkdir('{}'.format(todir))

    size=args.nproc
    with mp.Pool(processes=size) as pool:
        results = [pool.apply_async(getData,args=(i, size, todir,)) for i in range(size)]
        results_get = [r.get() for r in results]

    
    # Get overall error status
    es=0 # No error
    for e in results_get:
        if e==-1:
            es=-1 # Error in some processor

    print(es)
