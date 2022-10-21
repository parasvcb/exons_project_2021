const reactStringReplace = require("react-string-replace");

export const renderFasta = (transid, transexonids, exonHash) => {
//renderFasta(transid, transexonids, exonHash) {
    let fasta = "";
    let header = ">" + transid;
    let arr = [];
    arr.push(header);
    transexonids.split(",").map(ex => (fasta += exonHash[parseInt(ex)].aaseq));
    let uppfac = parseInt(fasta.length / 80);
    for (var i = 0; i <= uppfac; i++) {
      arr.push(fasta.slice(i * 80, (i + 1) * 80));
    }
    //console.log("fastaview,", arr);
    return (
      <span style={{ fontSize: 15, whiteSpace: "pre-line" }}>
        {arr.join("\n")}
      </span>
    );
}

export const transSisdomss = (looper, trid, property) => {
//transSisdomss(looper, trid, property) {
    //console.log("looper,", property);
    let tstate = true;
    if (looper.length < 2) {
      //console.log(looper.length === 0 ? "false" : looper[0][property]);
      return looper.length === 0 ? false : looper[0][property];
    } else {
      let patt = RegExp(`${trid}`, "g");
      for (var i = 0; i < looper.length; i++) {
        if (patt.exec(looper[i].list_trans_fk)) {
          //console.log(looper[i][property]);
          //console.log(property,looper[i][property])
          return looper[i][property];
        }
      }
    }
    if (tstate) {
      //console.log("false");
      return false;
    }
}

export const propConcate = (transfk, domset, singleexon, visprop) =>{
    //console.log("sing",singleexon)
    let transfkstr = transfk.toString();
    if (visprop=="domseq"){
    return transSisdomss(singleexon.exonsDom, transfk, visprop);
    }
    else if (visprop=="ssseq"){
      return transSisdomss(singleexon.exonsSS, transfk, visprop);
    }
    else {
      return transSisdomss(singleexon.exonsDis, transfk, visprop);
    }
}

export const decideText = (state, transfk, domset, singleexon) => {
    let transfkstr = transfk.toString();
    if (state.showAaseq) {
      return decorateText(
        singleexon.aaseq,
        singleexon.exId,
        [],
        "aaseq",
        []
      );
    }
    if (state.showSs) {
      let ss = transSisdomss(singleexon.exonsSS, transfk, "ssseq");
      return decorateText(ss, singleexon.exId, ["H", "E"], "ss", []);
    }
    if (state.showDis) {
      let dis = transSisdomss(singleexon.exonsDis, transfk, "disseq");
      return decorateText(dis, singleexon.exId, ["D"], "dis", []);
    }

    if (state.showDom) {
      let dom = transSisdomss(singleexon.exonsDom, transfk, "domseq");
      return decorateText(dom, singleexon.exId, [], "dom", domset);
    }
}

export const decorateText = (textassum, id, vals, keyst, domset) => {
    //let opaval = id[0] === "R" ? 0.25 : id.split(".")[2] === "A" ? 0.25 : id.split(".")[2] === "F" ? 0.5  : 0.7;
    //console.log("idexon", id, opaval);
    const spanBackgroundUTR ={
      backgroundColor: "white",
      color: "grey",
      width: "100px",
      border: "1px solid black",
      wordWrap: "break-word",
      marginRight: "10px"
    }
    const spanBackground = {
        backgroundColor: "#f8f9fa", //#efefef and linen
        fontWeight: "lighter",
        color: "#00303f",
        border: "1px solid black",
        marginRight: "10px"
      
    }
    let replacedText = textassum.length > 0 ? textassum : false;
    //console.log(id, replacedText, vals, keyst);
    if (replacedText) {
      //console.log("yes replaed text");
      if (keyst === "ss") {
        vals.map(iter => {
          let color = iter === "H" ? "purple" : "yellow";
          let patt = new RegExp(`(${iter}+)`, "g");
          replacedText = reactStringReplace(replacedText, patt, (match, i) => (
            //console.log("ss", match, i, match + i, id),
            <span key={match + i} style={{ background: color, color: color }}>
              {match}
            </span>
          ));
        });
      } else {
        if (keyst === "dom") {
          domset.map(iter => {
            let patt = new RegExp(`(${iter.code}+)`, "g");
            replacedText = reactStringReplace(
              replacedText,
              patt,
              (match, i) => (
                <span
                  key={match + i}
                  style={{
                    background: iter.color,
                    color: iter.color,
                    // border: "double 1px blue"
                  }}
                >
                  {match}
                </span>
              )
            );
          });
        } else {
          if (keyst === "dis") {
            let iter = "D";
            let patt = new RegExp(`(${iter}+)`, "g");
            replacedText = reactStringReplace(
              replacedText,
              patt,
              (match, i) => (
                //console.log("ss", match, i, match + i, id),
                <span
                  key={match + i}
                  style={{ background: "olive", color: "olive" }}
                >
                  {match}
                </span>
              )
            );
          }
        }
      }
    }
    if (replacedText) {
      //console.log("yesitis", replacedText, id);

      return (
        <span
          className="badge"
          style={spanBackground}
        >
          {replacedText}
        </span>
      );
    } else {
      return (
        <span
          className="badge"
          style={spanBackgroundUTR}
        >
          UTR
        </span>
      );
    }
}
