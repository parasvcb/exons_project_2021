import React, { Component } from "react";
import pic1 from "../pic_1.png";
import pic2 from "../pic_2.png";
import pic3 from "../pic_3.png";
import pic4 from "../pic_4.png";
import pic5 from "../pic_5.png";
import pic6 from "../pic_6.png";

/*
this tutorial will have brief description of the components and views of page, 
H2: tutorial section
     In Brief:
   database will render gene and its protein coding transctipts (NCBI sourced). Current functionalities of this will includes
       0) showing the nomenclature of the exon and depicting exon s aof varied nature and context in portein transcripts.
       a) showig exon beakup of the portein segments
       b) secondary structure, domain and disorder prediction view along with exon components of proteins. 
       c) showing exonic placeholder comparison of the protein coding transcripts.
   above functionalities will be represented and shown in the following views. 
      a) Exon Detailed view: with focus on the nomenclature per exon view,
      b) exon nightingale component 
      c) 
*/

//by default reference isoform will be renederd by the  
class Tutorial extends Component {
  
  render() {
    let header2 = "Tutorial Section"
    return (
      <div className="container">
        <div className="jumbotron">
          <br></br>
          <h2> {header2} </h2>
          <p>
            From Home page, gene can be queried using either the Gene ID (Entrez
            NCBI) or using its name.
          </p>
          <p>
            After selecting the Gene from the listed drop down list, following
            information will be rendered:
          </p>

          <ul>
            <li>
              The result page will be headed by the <b>Gene name</b> and{" "}
              <b>organism</b>. Information about the transcripts is given in the
              'TRANSCRIPT TABLE'. For each transcript, the number of exons
              belonging to different categories of 1<sup>st</sup> level
              classification (UTMRD) and 2<sup>nd</sup> level of classification
              (AGF) is also given in the table.
              <br />
              <br />
              <img
                src={pic1}
                style={{ height: "75%", width: "75%", objectFit: "contain" }}
              />
            </li>
            <br />
            <li>
              'EXON TABLE' contains the detailed description of every exon in
              all the transcripts along with their descriptive nomenclature
              (Exon Id), Raw and Coding genomic coordinates (RGSC,RGEC - Raw
              Genomic start/end coordinates; CGSC,CGEC- Coding genomic start/end
              coordinates), length of coded amino acids and weighted exon
              frequency (WEF).
              <br />
              <br />
              <img
                src={pic2}
                style={{ height: "75%", width: "75%", objectFit: "contain" }}
              />
            </li>
            <br />
            <li>
              The color coded, Pfam annotated domains can be seen in the 'DOMAIN
              TABLE'. After this table, all the transcripts are listed and their
              details can be fetched by clicking on 'SHOW' button.
            </li>
            <br />
            <img
              src={pic3}
              style={{ height: "75%", width: "75%", objectFit: "contain" }}
            />
            <br />
            <br />
            <li>
              After clicking on 'SHOW' button against the transcript, exon
              details can be seen. For non-coding exons, content is displayed by
              #, whereas for coding exons it is single letter amino acid code.
            </li>
            <br />
            <img
              src={pic4}
              style={{ height: "75%", width: "75%", objectFit: "contain" }}
            />
            <br />
            <br />
            <li>
              The secondary structure of the exons can be highlighted on the
              exons by clicking on 'Show SS', yellow is 'Beta sheet region',
              purple is 'Alpha helix region' and '-' coiled/unordered region.
            </li>

            <br />
            <img
              src={pic5}
              style={{ height: "75%", width: "75%", objectFit: "contain" }}
            />
            <br />
            <br />
            <li>
              The Pfam domains of the exons can be visualised by clicking on
              'Show Dom', Uncolored region not colored or in '-' is not coding
              any domains.
            </li>
            <br />
            <img
              src={pic6}
              style={{ height: "75%", width: "75%", objectFit: "contain" }}
            />
            <br />
          </ul>
        </div>
      </div>
    );
  }
}
const headerStyle = {
  background: "Seagreen",
  fontWeight: 500,
  color: "Honeydew",
  padding: "1rem",
  fontFamily: "font-family: Arial, Helvetica, sans-serif"
};
const boldStyle = {
  color: "white",
  fontWeight: "thick"
};

export default Tutorial;
