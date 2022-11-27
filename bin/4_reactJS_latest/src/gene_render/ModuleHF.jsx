//exonalerts(exons) {
// change above too follwing to make it work
export const MDexonalerts = (exons) =>{
    let ob = {};
    exons.forEach(element => {
      let messageComp1 = "";
      if (element.exId[0] !== "R") {
        const id = element.exId.split(".");
        let message1h = id[0] + "." + id[1];
        let message1 = "";
        switch (id[0]) {
          case "U":
            message1 +=
              "This exon is always part of UTR in all the isoforms(ISFs) and always non-coding ";
            break;
          case "T":
            message1 +=
              "This exon is always part of CDS in all the isoforms(ISFs) ";
            break;
          case "M":
            message1 +=
              "This exon is always part of CDS in all the isoforms(ISFs) but has only 1nt contribution and hence no Amino acid ";
            break;
          case "D":
            message1 +=
              "This exon is either part of CDS or UTR in different isoforms(ISFs) of this gene ";
            break;
          default:
            message1 += " ";
        }
        if (
          ["D-1", "D-2", "D0", "T-1", "T-2", "T0"].indexOf(id[0] + id[1]) > -1
        ) {
          switch (id[0] + id[1]) {
            case "D-1":
              message1 +=
                "and it is non-coding due to PTC in upstream region for this ISF;";
              break;
            case "D-2":
              message1 += "and it is non-coding in this ISF;";
              break;
            case "D0":
              message1 += "and it has only 1nt share to CDS in this ISF;";
              break;
            case "T-1":
              message1 +=
                "but is not coding for aa in this ISF due PTC in upstream region;";
              break;
            case "T0":
              message1 += "but it has only 1nt share to CDS for this ISF;";
              break;
          }
        } else {
          switch (id[1]) {
            case "-2":
              message1 += "";
              break;
            case "-1":
              message1 += "";
              break;
            case "0":
              message1 += "";
              break;
            case "1":
              message1 += " and is coding amino acids (Reference) in this ISF,";
              break;
            default:
              message1 +=
                " and is coding amino acids different than that of reference exon (Variant:" +
                id[1] +
                ") in this ISF,";
          }
        }

        let message2h = id[2] + "." + id[3];
        let message2 = " further this exon is ranked ";

        if (id[3].length > 1 && id[3][id[3].length - 2] === "1") {
          message2 += id[3] + "th";
        } else {
          switch (id[3].slice(-1)) {
            case "1":
              message2 += id[3] + "st";
              break;
            case "2":
              message2 += id[3] + "nd";
              break;
            case "1":
              message2 += id[3] + "rd";
              break;
            default:
              message2 += id[3] + "th";
              break;
          }
        }

        message2 += " in gene sequence and ";
        switch (id[2]) {
          case "G":
            message2 += " is constitutively present in all the ISFs,";
            break;
          case "A":
            message2 +=
              " is alternatively present in ISFs with WEF:" + element.wef + ",";
            break;
          case "F":
            message2 +=
              " is present in all the ISFs with certain alternative splice site variations,";
            break;
          default:
            message2 += "";
        }

        let message3h = id[4] + "." + id[5];
        let message3 = "";
        switch (id[4]) {
          case "0":
            message3 +=
              " no alternative splice site is choosen for this exon in this ISF";
            break;
          case "n":
            message3 +=
              " alternative 5' splice site has been choosen for this exon different from reference exon ";
            break;
          case "c":
            message3 +=
              " alternative 3' splice site has been choosen for this exon different from reference exon ";
            break;
          case "b":
            message3 +=
              " alternative 5' splice site and alternative 3' splice site has been choosen which are different from reference exon ";
            break;
          default:
            message3 += "";
        }
        if (["n", "c", "b"].indexOf(id[4]) > -1) {
          message3 += element.parent;
          message3 += " and that is " + id[5] + "such occurence";
        }
        messageComp1 = (
          <div>
            <p style={{ display: "inline", fontWeight: "bold" }}>
              <span style={{ color: "crimson" }}>{id[0] + "." + id[1]}.</span>
              <span style={{ color: "DarkGreen" }}>{id[2] + "." + id[3]}.</span>
              <span style={{ color: "DarkOrange" }}>
                {id[4] + "." + id[5]}.
              </span>
            </p>
            <p style={{ fontSize: 15 }}>
              <span style={{ color: "crimson", fontWeight: "bold" }}>
                ({message1h}){" "}
              </span>
              <span>{message1}</span>
              <span style={{ color: "DarkGreen", fontWeight: "bold" }}>
                ({message2h}){" "}
              </span>
              <span>{message2}</span>

              <span style={{ color: "DarkOrange", fontWeight: "bold" }}>
                ({message3h}){" "}
              </span>
              <span>{message3}</span>
            </p>
          </div>
        );
      } else {
        let id = element.exId.split(":");
        let message1h = id[0] + ":" + id[1];
        let message2h = id[2];
        let message3h = id[3];
        let message4h = id[4];
        let message1 = "This is a intron retention case which";
        switch (id[1]) {
          case "-2":
            message1 += " is non-coding in this transcript";
            break;
          case "1":
            message1 += " is coding in this transcript ";
            break;
          case "0":
            message1 += " has only 1nt contribution and hence no amino acid ";
            break;
          case "-1":
            message1 +=
              " is not coding for aa in this transcript due PTC in upstream region ";
            break;
          default:
            message1 += " ";
        }
        let message2 = " and spans from exon ";
        let message3 = " and is the ";
        let message4 = " and ending at exon ";
        messageComp1 = (
          <div>
            <strong>
              <p style={{ display: "inline", fontWeight: "bold" }}>
                <span style={{ color: "crimson" }}>{id[0] + ":" + id[1]}:</span>
                <span style={{ color: "DarkGreen" }}>{id[2]}:</span>
                <span style={{ color: "DarkOrange" }}>{id[3]}:</span>
                <span style={{ color: "SteelBlue" }}>{id[4]}</span>
              </p>
            </strong>
            <p>
              <span style={{ color: "crimson", fontWeight: "bold" }}>
                {message1h}{" "}
              </span>
              <span>{message1}</span>
              <span>{message2}</span>
              <span style={{ color: "DarkGreen", fontWeight: "bold" }}>
                {id[2]}
              </span>
              <span>{message3}</span>
              <span style={{ color: "DarkOrange", fontWeight: "bold" }}>
                {id[3]}
              </span>
              <span> instance starting from this exon</span>
              <span>{message4}</span>
              <span style={{ color: "SteelBlue", fontWeight: "bold" }}>
                {id[4]}
              </span>
            </p>{" "}
          </div>
        );
      }
      ob[element.exId] = messageComp1;
    });
    return ob;
}


export const MDlastexon = (exons) => {
  //lastexon(exons) {
    exons.sort(function(a, b) {
      return a.rawst - b.rawst;
    });
    const lastele = exons[exons.length - 1];
    if (lastele != "R") {
      return parseInt(lastele.exId.split(".")[3]);
    } else {
      return parseInt(lastele.exId.split(":")[2].split(".")[3]);
    }
    // returns the numeric placeholder of the last exon
}

export const MDexoncodes = (exons) => {
    let finalob = {};
    exons.forEach(element => {
      if (element.exId[0] !== "R" && element.exId.split(".")[4] === "0") {
        finalob[parseInt(element.exId.split(".")[3])] = element.exId[0];
      }
    });
    //console.log("EXONCODES", finalob);
    return finalob;
    //1: "T",2: "U",3: "T",4: "T" output
    //excluding retenetion cases and that of the 0 status exons[4th element], stiore their numeric cplaceholders as key and global coding state as value.
}

// export const MDaaSeqPerTrans =(ob,exonHash) => {
//     let aaVista=''
//     ob.trans.map(trans => {
//       exonOb=[] // this should have individual exon entry with name of exon, 
//       let aaVista=''
//       let sseqVista=''
//       let domseqVista=''
//       let disseqVista=''
//       let exons={}
//       trans.exonsIds.split(",").map(ex => {
//         aaVista += exonHash[parseInt(ex)].aaseq
//         sseqVista += propConcate(trans.id,domset,exonHash[parseInt(ex)],"ssseq")
//         domseqVista += propConcate(trans.id,domset,exonHash[parseInt(ex)],"domseq")
//         disseqVista += propConcate(trans.id,domset,exonHash[parseInt(ex)],"disseq")
//         //aavista should be sequqnce IDentifier for this trans, keep them added in cases viewing by coponent is needed
//         // exon attribute will be ID, locations:[{fragments:[start:, end:,shape:]}]
//         //domains will be 
//         //use other to find regex
//       })});
// }


function union(setA, setB) {
  const _union = new Set(setA);
  for (const elem of setB) {
    _union.add(elem);
  }
  return _union;
}

function intersection(setA, setB) {
  const _intersection = new Set();
  for (const elem of setB) {
    if (setA.has(elem)) {
      _intersection.add(elem);
    }
  }
  return _intersection;
}


const dictionaryHexColor = {
  1.00  : "FF",
  0.99  : "FC",
  0.98  : "FA",
  0.97  : "F7",
  0.96  : "F5",
  0.95  : "F2",
  0.94  : "F0",
  0.93  : "ED",
  0.92  : "EB",
  0.91  : "E8",
  0.90  : "E6",
  0.89  : "E3",
  0.88  : "E0",
  0.87  : "DE",
  0.86  : "DB",
  0.85  : "D9",
  0.84  : "D6",
  0.83  : "D4",
  0.82  : "D1",
  0.81  : "CF",
  0.80  : "CC",
  0.79  : "C9",
  0.78  : "C7",
  0.77  : "C4",
  0.76  : "C2",
  0.75  : "BF",
  0.74  : "BD",
  0.73  : "BA",
  0.72  : "B8",
  0.71  : "B5",
  0.70  : "B3",
  0.69  : "B0",
  0.68  : "AD",
  0.67  : "AB",
  0.66  : "A8",
  0.65  : "A6",
  0.64  : "A3",
  0.63  : "A1",
  0.62  : "9E",
  0.61  : "9C",
  0.60  : "99",
  0.59  : "96",
  0.58  : "94",
  0.57  : "91",
  0.56  : "8F",
  0.55  : "8C",
  0.54  : "8A",
  0.53  : "87",
  0.52  : "85",
  0.51  : "82",
  0.50  : "80",
  0.49  : "7D",
  0.48  : "7A",
  0.47  : "78",
  0.46  : "75",
  0.45  : "73",
  0.44  : "70",
  0.43  : "6E",
  0.42  : "6B",
  0.41  : "69",
  0.40  : "66",
  0.39  : "63",
  0.38  : "61",
  0.37  : "5E",
  0.36  : "5C",
  0.35  : "59",
  0.34  : "57",
  0.33  : "54",
  0.32  : "52",
  0.31  : "4F",
  0.30  : "4D",
  0.29  : "4A",
  0.28  : "47",
  0.27  : "45",
  0.26  : "42",
  0.25  : "40",
  0.24  : "3D",
  0.23  : "3B",
  0.22  : "38",
  0.21  : "36",
  0.20  : "33",
  0.19  : "30",
  0.18  : "2E",
  0.17  : "2B",
  0.16  : "29",
  0.15  : "26",
  0.14  : "24",
  0.13  : "21",
  0.12  : "1F",
  0.11  : "1C",
  0.10  : "1A",
  0.09  : "17",
  0.08  : "14",
  0.07  : "12",
  0.06  : "0F",
  0.05  : "0D",
  0.04  : "0A",
  0.03  : "08",
  0.02  : "05",
  0.01  : "03",
  0.00  : "00"
}
// https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/Set

export const MDaddNightingaleReadyJSON = (ob, exonAlert, exonHash, exonWif) => {  
    // shapes can be roundRectangle, discontinuosStart, discontinuos, discontinuosEnd, rectangle, line, bridge, diamond, chevron, catFace, triangle, wave, hexagon, pentagon, circle, arror, doubleBar, helix, strand
    ////////////////////////////////////////////////////////////////////////////////
    
    /*this () will take gene ob JSON, exonAlert, defined above and called from  _main, and exonHash having id to exon Object. 
    We need to add new JSON (for nightingale viz) in nightingale format, and modify (add) them to individual transcrip objects in ob.trans[]. We will be adding exons (interpro domains format, f:'interpro.js', coods aa), 
    secondary structures and disorder (ss format,f:'interpro-secondary-structure.json') format, and domains in domains format 
    */

    /*  logic: iterate the trans ob (list of transcript objects), define list of nightingale components and their objects (4 are defined below). 
    
    iterate then the exons list, and for every exon list, call exon ob using exon Hash. (record Sequence also)

                    iterate exons, and push them to nihtingale component (only that have aa>0 aa), other to nightingale description datatable
    */
    // console.log(ob)
    // console.log("OB",ob.trans)
    // console.log(typeof(ob.trans))
    // console.log(Array.isArray(ob['trans']))
    // console.log("PRINT OB", ob)
    // console.log("PRINT exonHash", exonHash)
    ob.trans.map(trans => {
              //console.log("inside")
              let exonNightingale=[] // this should have individual exon entry with name of exon, 
              // will fetch {} objects with accession, locations[{...}], and color
              let ssNightingale=[
                {
                  "accession":trans.tId,
                  "locations":[
                    {
                      "fragments":[]
                    }
                  ]
                }
              ]
              let disorderNightingale=[
                {
                  "accession":trans.tId,
                  "locations":[
                    {
                      "fragments":[]
                    }
                  ]
                }
              ]
              // list of 
              let domainsNightingale=[]
              // akin to exons
              let dataTableNightingale=[]
              // {begin,end,type,description}
              let aaVista=''
              let domStringNightingale = "overlapping"
              let spanDomainSet = new Set ()
              // console.log("going In")
              trans.exonsIds.split(",").map(ex => {
                // console.log("ex",ex)
                aaVista += exonHash[parseInt(ex,10)].aaseq
               })
              // console.log("going Out")
              trans.sequence=aaVista
              //sequence is recorded
              /////////////////////////////////////////////////////////////////////////
              //     [ {  "accession": "chain_A", "locations": [
              //      {
              //        "fragments": [
              //          {
              //            "start": 45, "end": 100, "shape": "helix", "fill": "transparent", "color": "#f06"
              //          },{... sheet/helix and more} ] } ] } ]
              //"exonsRegion": "U.-2.A.4.n.1$_U.-2.A.6.c.1$_T.1.A.6.n.1$1,240,
              let countExon=1
              trans.exonsRegion.split('_').map(ex => {
                          let elem=ex.split('$')
                          //elem[0]=id, 1 is 1,240,shape
                          let coods=elem[1].split(',')
                          let start=parseInt(coods[0])
                          let end=parseInt(coods[1])
                          let shapeInt=parseInt(coods[2])
                          let shape="roundRectangle"
                          if (shapeInt==1){
                            shape="discontinuosStart"
                          } else if (shapeInt==2){
                            shape="discontinuosEnd"
                          } else if (shapeInt==3){
                            shape="discontinuos"
                          } else {
                            shape="rectangle"
                          }
                          // pending: else above will nenevr render because said condition wont ever turn true
                          if (start!=end){
                          //only push if they arent same
                            // let colExon = (countExon %2 == 0) ? "grey" : "black" // "#CCCCCC" was default 
                            let colExon = (countExon %2 == 0) ? "#1A0315" : "#88D317" // "#deepplum, electric lime" was default 
                            let wef = exonWif[elem[0]]
                            let val2 = Math.round(parseFloat(wef)*1e2)/1e2;
                            // // if you need 3 digits, replace 1e2 with 1e3 etc.
                             
                            // console.log (typeof(val2), elem[0], wef, val2, dictionaryHexColor[val2])
                            colExon = colExon + dictionaryHexColor[val2]
                            let exonNight = {"accession":elem[0], "color": colExon, 'locations':[{
                                'fragments':[{'start':start, 'end': end, 'shape':shape}]}]}
                            countExon += 1
                            // console.log("INTERNAL MOD", coods, countExon, colExon)
                            // interpro.js
                          exonNightingale.push(exonNight)
                          }
                        let exonDat={'start':start, 'end': end,'type':'exons','description':elem[0], 'tooltip':exonAlert[elem[0]]}
                        if (start==end) {
                          exonDat['tooltip']=<><><h4>This exon is part of UTR, and has gene coordinates {coods[3]}-{coods[4]}, this is not rendered above.</h4></><>exonAlert[elem[0]]</></>
                        }
                        dataTableNightingale.push(exonDat)
              })
              console.log("1st pass DONE")
              trans.secondaryStructure.split('-').map(hesTypes => {
                        //"H:9,15_37,47-E:221,225_351,357_-S:1,4
                        //H: E: S: tags will be there always but there entity after : may be missing so be a little careful.
                        // hstypes will have string of secondarys trtcutres,
                        let indSS=hesTypes.split(':')
                        // first element will be H and second will be spans list
                        if (hesTypes.length>2) {
                          //"H:", this length will be still 2, span should have length of 4 and greater than that.
                          // only exon coordinates, UTR entity, shall we have genomic coods or not

                          let eleToPush={"start":'',"end":'',"shape":"helix","fill":"transparent","color":"#f06"}
                          let datString="Helical region prediction"

                          if (indSS[0]=="E") {
                            eleToPush={"start":'',"end":'',"shape":"strand","color":"#fc0"}
                            datString="Beta Strand prediction"
                          } else if (indSS[0]=="D") {
                            eleToPush={"start":'',"end":'',"shape":"rectangle","fill":"#565656","color":"#565656"}
                            datString="Disorder region prediction"
                          }

                          indSS[1].split('_').map(spans => {
                            let spanEle=spans.split(',')
                            let start=parseInt(spanEle[0])
                            let end=parseInt(spanEle[1])
                            let ssNight= { ...eleToPush}
                            ssNight.start=start
                            ssNight.end=end
                            let ssDat=''
                            if (indSS[0]=="D") {
                              disorderNightingale.push(ssNight)
                              ssDat={'start':start, 'end': end,'type':'Disorder region','description':datString, 'tooltip':''}
                            } else {
                              // ssElements
                              ssNightingale.push(ssNight)
                              ssDat={'start':start, 'end': end,'type':'Secondary structure','description':datString, 'tooltip':''}  
                            }
                            dataTableNightingale.push(ssDat)
                          })}})
              
              trans.domains.split('$').map(domains => {
                        //z-alpha,PF02295.17,aqua:1,63$PF00035.26,lawngreen:207,271_318,382_432,496
                        //NoDomainEntryPredicted,0,black:0,0_
                        let comp=domains.split(':')
                        let spans=comp[1]
                        let elem=comp[0].split(',')
                        let name=elem[0]
                        let pfamId=elem[1]
                        let color=elem[2]
                        spans.split('_').map(span => {
                              let spanEle=span.split(',')
                              let start=parseInt(spanEle[0])
                              let end=parseInt(spanEle[1])

                              if (domStringNightingale === "overlapping") {
                                let smallset = new Set()
                                for (let ti=start; ti <= end; ti++) {
                                  smallset.add(ti)
                                }
                                const inter = intersection (smallset, spanDomainSet)

                                if (inter.size >0) {
                                  domStringNightingale = "non-overlapping"
                                }
                                else {
                                  spanDomainSet = union (smallset, spanDomainSet)
                                }
                              }
                              
                              // spanDomainSet

                              let domNight = {
                                "accession":pfamId, 
                                "color": color, 
                                'locations':[
                                  {
                                    'fragments':
                                      [
                                        {
                                          'start':start, 
                                          'end': end, 
                                          'shape':'roundRectangle'
                                        }
                                      ]
                                  }
                                ]}
                              domainsNightingale.push(domNight)
                              let domDat={
                                'start': start, 
                                'end': end, 
                                'type':'domains', 
                                'description':pfamId, 
                                'tooltip':<>
                                          <h5 className="text-secondary">
                                            <span>Domain entry</span>&nbsp; 
                                            <span>
                                            <a className="text-decoration-none text-primary"
                                             target="_blank" 
                                            //  href={"http://pfam-legacy.xfam.org/family/"+pfamId}>
                                             href={"https://www.ebi.ac.uk/interpro/entry/pfam/"+pfamId.split('.')[0]}>
                                            {pfamId}
                                      
                                            </a> 
                                            </span>&nbsp;
                                            <span>{name}</span>
                                          </h5></>
                              }
                              dataTableNightingale.push(domDat)
                              console.log (trans.tId)
                              console.log (domNight)
                              console.log (domDat)
                              console.log ("&&&&&")
                          })
                        })

                          //list1.sort((a, b) => (a.start > b.start) ? 1 : -1)
                          dataTableNightingale.sort((a, b) => (a.start > b.start) ? 1 : -1)
                          trans['nightSs']=ssNightingale
                          trans['nightDis']=disorderNightingale
                          trans['nightDom']=domainsNightingale
                          trans['nightExons']=exonNightingale
                          trans['nighthTable']=dataTableNightingale
                          trans['domnightLayout'] = domStringNightingale
                          //console.log("trans--<",trans)
    })
    // console.log("Returned ??")
    return ob;
};

export const MDgeneCard = (gene) => {
  //geneCard(gene) {
    let cod = 0;
    let non = 0;
    gene.exonsGenes.forEach(element => {
      if (element.aaseq.length > 0) {
        cod = cod + 1;
      } else {
        non = non + 1;
      }
    });
    return (
      // changing them from cards to gutters because they should be of equal height, width i am changing
      
        <div className="row gx-5 mt-1 gy-5">
          {/* "p-2 w-75 bg-light" */}
          <div className="col-md-9 m-0">
          <div className="col-md-12 m-0 p-3" style={cardStyle1}>
              <a
                className="text-decoration-none"
                title="To NCBI Gene Page"
                placement="auto"
                rel="noreferrer"
                target="_blank"
                href={"https://www.ncbi.nlm.nih.gov/gene/" + gene.entrezid}
              >
                <h4 className="m-0 p-0 text-black">
                  Gene: {gene.name}
                </h4>
              </a>
              <br></br>
              <a
                  target="_blank"
                  rel="noreferrer"
                  title="To NCBI Gene Page"
                  placement="auto"
                  className="text-decoration-none"
                  href={"https://www.ncbi.nlm.nih.gov/gene/" + gene.entrezid}
              >
                  <p className="m-0 p-0 text-black ">NCBI ID: {gene.entrezid}</p>
              </a>
              
              <a
                  target="_blank"
                  rel="noreferrer"
                  title="NCBI taxonomy link"
                  placement="auto"
                  className="text-decoration-none"
                  href={
                    "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" +
                    gene.txid
                  }
                >
                <p className="m-0 pb-2 text-black ">Organism: {gene.organism}</p>
              </a>
            </div>
            </div>
        <div className="col-md-3 m-0">
        <div className="col-md-12 m-0 p-3" style={cardStyle2}>
          <div>
            <h4 >Constituents</h4>
              <div >
                <p className="m-0 p-0">Transcripts: {gene.trans.length}</p>
                <p className="m-0 p-0">Coding Exons: {cod}</p>
                <p className="m-0 p-0">NonCoding Exons: {non}</p>
              </div>
          </div>
        </div>
        </div>
        </div>
    );
  }
  const cardStyle1 = {
    background: "#e9e9e9",//DCAE1D, 494c6b
    // background:"linear-gradient(to right, #d5d5d5, #e9e9e9)",
    color: "black",
    borderRadius: "20px",
    // padding: "1px",
    // boxShadow: "10px 10px 5px lightblue",
    // margin:"0px 1px 1px 0px"
  };

  const cardStyle1_1 = { 
      backgroundColor: "transparent",
      backgroundImage: "linear-gradient(to bottom, rgba(30, 87, 153, 0.2) 0%, rgba(125, 185, 232, 0) 100%)",
      backgroundRepeat: "repeat",
      borderRadius:"10px",
      padding: "5px"
  
  }
  const cardStyle2 = {
    background: "#e9e9e9",//7A9D96, 494c6b
    // background:"linear-gradient(to right, #e9e9e9, #d5d5d5)",
    color: "black",
    borderRadius: "25px",
    // padding: "1px"
    // paddingBottom: "10" //this wont get executed bro
  
  };