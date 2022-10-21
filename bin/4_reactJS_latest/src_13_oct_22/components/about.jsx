import React, { Component } from "react";
import about from "../pic_about.png";
class About extends Component {
  render() {
    return (
      <div className="container overflow-scroll">
          <br></br>
          <h2>Overview</h2>
          <p>
            NEXTRaP-db is a repository which gives detailed information of the
            structural overview of the transcripts of genes of 5 model
            organisms;{" "}
            <i>
              Homo sapiens, Mus musculus, Caenorhabditis elegans, Danio rerio
            </i>{" "}
            and <i>Drosophila melanogaster</i> taken from NCBI. Search can be
            done using, either the gene name or the EntrezID. The user can trace
            down the differences amongst different transcripts of the gene, both
            at genomic as well as proteomic level, using the exonic
            representation of the protein coded by each transcript. For every
            gene, following information is there in the form of tables:
          </p>
          <ul>
            <li>
              <b>TRANSCRIPT TABLE</b>: This table consists of information about
              each transcript regarding the domains, number and type of exons it
              contains. It also depicts which transcript was chosen as PI
              (Principal Isoform). The table is sorted according to the coding
              length of transcripts.{" "}
            </li>
            <li>
              <b>EXON TABLE:</b> It gives the detailed information about all the
              exons that are contained in the aforementioned transcripts, sorted
              according to the Raw Genomic Coordinates:
            </li>

            <ul>
              <li>
                {" "}
                <u>Exon Id</u>: Describes the name given to the exon according
                to the current nomenclature.{" "}
              </li>
              <li>
                <u>Raw Genomic Start Coordinates (RGSC)</u>: It refers to the
                starting point of the actual coordinates of the exon in the
                genome.
              </li>
              <li>
                <u>Raw Genomic End Coordinates (RGEC)</u>: It refers to the
                ending point of the actual coordinates of the exon in the
                genome.
              </li>
              <li>
                <u>Start Coding Genomic Coordinates (SCSC)</u>: It refers to the
                starting point of the coding coordinates (obtained by mapping
                the amino acid sequence to the coding sequence by making the 3n
                pairs) of the exon.
              </li>
              <li>
                <u>End Coding Genomic Coordinates (ECEC)</u>: It refers to the
                ending point of the coding coordinates (obtained by mapping the
                amino acid sequence to the coding sequence by making the 3n
                pairs) of the exon.
              </li>
              <li>
                <u>Length (aa):</u> Length of the amino acids coded by that
                exon.
              </li>
              <li>
                <u>Weighted Inclusion Frequency (WEF):</u> WEF of each exon is
                calculated by fractioning the number of transcripts that
                includes this exon by total transcripts.
              </li>
            </ul>
          </ul>

          <li>
            <b>DOMAIN TABLE</b>: It gives the color coded information of the
            Pfam domains (if annotated) present in the transcripts.
          </li>
          <li>
            <b>TRANSCRIPT INFORMATION:</b> The last section comprises of the{" "}
            <b>Exon Id</b> and <b>Exon Details</b>
            of every transcript. On clicking the <b>SHOW</b> button against
            every transcript Id, user can view its constituent exons arranged
            vertically in the form of boxes, which encapsulate the amino acid
            sequence they code for, with their assigned names (according to the
            nomenclature). The color of the boxes containing Exon Ids depict the
            type of exon, as defined in the first level of nomenclature.
            <br />
            <div className="container" align="middle">
              <img
                src={about}
                style={{ height: "40%", width: "40%", objectFit: "contain" }}
              />
            </div>
            <br />
            By clicking on the <b> SHOW SS</b> button, the exons will display
            their color coded secondary structure information where purple
            corresponds to <b>alpha helix</b>, yellow corresponds to{" "}
            <b>beta sheet</b> and hyphen corresponds to <b>coil/unstructured</b>{" "}
            region. The <b>SHOW Dom </b>button will reveal the local presence of
            transcript domains per exon, represented through colors as listed in
            the DOMAIN TABLE, in all the exons of the transcript. The amino acid
            information of the exons can be fetched back by clicking on the{" "}
            <b>Show AA </b>
            button.
          </li>
        </div>
    );
  }
}
export default About;
