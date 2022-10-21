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

export const MDaddNightingaleReadyJSON = (ob, exonAlert, exonHash) => {
    //shapes can be roundRectangle, discontinuosStart, discontinuos, discontinuosEnd, rectangle, line, bridge, diamond, chevron, catFace, triangle, wave, hexagon, pentagon, circle, arror, doubleBar, helix, strand

    ////////////////////////////////////////////////////////////////////////////////
    // this () will take gene ob JSON, exonAlert, defined above and called from   //
    // _main, and exonHash having id to exon Object. We need to add new JSON in   //
    // nightingale format, and modify (add) them to indivdual transcrip objects in//
    // ob.trans[]. We will be adding exons (interpro domains format, f:'interpro. // 
    // js', coods aa), secondary structures and disorder (ss format,f:            //
    // 'interpro-secondary-structure.json') format, and domains in domains format //
    /*  logic: iterate the trans ob (list of transcript objects), define list of nightingale components and their objects (4 are defined below). 
                    iterate then the exons list, and for every exon list, call exon ob using exon Hash. (record Sequence also)

                    iterate exons, and push them to nihtingale component (only that have aa>0 aa), other to nightingale description datatable
    */
    // console.log(ob)
    // console.log("OB",ob.trans)
    // console.log(typeof(ob.trans))
    // console.log(Array.isArray(ob['trans']))
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
      // list of 
      let domainsNightingale=[]
      // akin to exons
      let dataTableNightingale=[]
      // {begin,end,type,description}
      let aaVista=''
      trans.exonsIds.split(",").map(ex => {aaVista += exonHash[parseInt(ex)].aaseq })
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
            shape="roundRectangle"
          }
          // pending: else above will nenevr render because said condition wont ever turn true
          if (start!=end){
          //only push if they arent same
            let exonNight = {"accession":elem[0], "color": "#CCCCCC", 'locations':[{
                'fragments':[{'start':start, 'end': end, 'shape':shape}]}]}
            // interpro.js
          exonNightingale.push(exonNight)
          }
        let exonDat={'begin':start, 'end': end,'type':'exons','description':exonAlert[elem[0]]}
        if (start==end) {
          exonDat['description']=<><><h2>This exon is part of UTR, and has gene coordinates{coods[3]}{coods[4]}, this is not rendered above.</h2></><>exonAlert[elem[0]]</></>
        }
        dataTableNightingale.push(exonDat)
      })
      
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

          if (indSS[0]=="E"){
            eleToPush={"start":'',"end":'',"shape":"strand","color":"#fc0"}
            datString="Beta Strand prediction"
          } else if (indSS[0]=="S") {
            eleToPush={"start":'',"end":'',"shape":"roundRectangle","fill":"transparent","color":"#000000"}
            datString="Disorder region prediction"
          }
          indSS[1].split('_').map(spans => {
            let spanEle=spans.split(',')
            let start=parseInt(spanEle[0])
            let end=parseInt(spanEle[1])
            let ssNight= { ...eleToPush}
            ssNight.start=start
            ssNight.end=end
            ssNightingale.push(ssNight)
            let ssDat={'begin':start, 'end': end,'type':'Secondary Structure','description':datString}
            dataTableNightingale.push(ssDat)
          })}})
      
      trans.domains.split('$').map(domains => {
        //z-alpha,PF02295.17,aqua:1,63$PF00035.26,lawngreen:207,271_318,382_432,496
        let comp=domains.split(':')
        let spans=comp[1]
        let elem=comp[0].split(',')
        let name=elem[0]
        let pfamId=elem[1]
        let color=elem[2]
        spans.split('_').map(span => {
          let spanEle=spans.split(',')
          let start=parseInt(spanEle[0])
          let end=parseInt(spanEle[1])
          let domNight = {"accession":pfamId, "color": color, 'locations':[{
            'fragments':[{'start':start, 'end': end, 'shape':'roundRectangle'}]}]}
          
          domainsNightingale.push(domNight)
          let domDat={'begin':start, 'end': end,'type':'domains','description':<><h2>Domain entry {pfamId}({name})</h2></>}
          dataTableNightingale.push(domDat)
          })})
          trans['nightSs']=ssNightingale
          trans['nightDom']=domainsNightingale
          trans['nightExons']=exonNightingale
          trans['nighthTable']=dataTableNightingale
          //console.log("trans--<",trans)
    })
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
      <div className="d-flex">
        <div className="p-2 w-75 bg-light">
          <div className="card" style={cardStyle1}>
            <a
              target="_blank"
              href={"https://www.ncbi.nlm.nih.gov/gene/" + gene.entrezid}
            >
              <h4 className="card-title m-0 p-0 text-white">
                Gene: {gene.name}
              </h4>
            </a>
            <div className="card-body p-2 mb-2 text-white">
              <a
                target="_blank"
                href={"https://www.ncbi.nlm.nih.gov/gene/" + gene.entrezid}
              >
                <p className="m-0 p-0 text-white">NCBI ID: {gene.entrezid}</p>
              </a>
              <a
                target="_blank"
                href={
                  "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" +
                  gene.txid
                }
              >
                <p className="m-0 pb-2 text-white">Organism: {gene.organism}</p>
              </a>
            </div>
          </div>
        </div>
        <div className="p-2 w-25 bg-light  ">
          <div className="card" style={cardStyle2}>
            <h4 className="card-title m-0 p-0">Constituents</h4>
            <div className="card-body p-1 m-0">
              <p className="m-0 p-0">Transcripts: {gene.trans.length}</p>
              <p className="m-0 p-0">Coding Exons: {cod}</p>
              <p className="m-0 p-0">NonCoding Exons: {non}</p>
            </div>
          </div>
        </div>
      </div>
    );
  }
  const cardStyle1 = {
    background: "steelblue",
  
    paddingBottom: "10",
    color: "white"
  };
  const cardStyle2 = {
    background: "dodgerblue",
    color: "white",
  
    paddingBottom: "10"
  };