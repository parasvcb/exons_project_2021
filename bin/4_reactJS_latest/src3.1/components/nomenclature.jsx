import React, { Component } from "react";

class Nomenclature extends Component {
  render() {
    return (
      <div className="container">
        <div className="jumbotron">
          <h2>Exon nomenclature</h2>
          <h5>
            The exon nomenclature consists of 6 block code separated by (.) and
            is divided into following three-level alphanumeric categories
            detailing various attributes of exons.
          </h5>
          <br />

          <h3 align="center">
            <span style={{ color: "SeaGreen" }}>[UTMDR].[(-2)-n].</span>
            <span style={{ color: "MediumVioletRed" }}>[GAF].[1-n].</span>
            <span style={{ color: "Teal" }}>[ncb0].[0-n]</span>
          </h3>
          <br />
          <ul>
            <li>
              <h5 style={{ color: "SeaGreen" }}>
                Level 1 -> To describe the correlation of exons to the
                translated product
              </h5>
              <br /> <b>U:</b> Exons, which remain part of Untranslated region
              UTR in every transcript; numeric code is -2
              <br />
              <b>T</b>: Exons, which are always part of Coding sequence CDS
              region; numeric code is -1 to n
              <ul>
                <li>
                  {" "}
                  -1 depicts, premature stop codon in upstream exon and hence
                  this cannot code for amino acid
                </li>
                <li>
                  {" "}
                  0, means it is part of CDS but CGC span is only 1 nucleotide
                </li>
                <li>
                  1 means it is contributing amino acids in transcripts, and
                </li>
                <li>
                  {" "}
                  2 or more depicts, this exon has varying/different amino acids
                  contribution from transcript to transcript (resulting from FSE
                  in upstream exon or ATI/ATT in this exon)
                </li>
              </ul>
              <b>M:</b> Exons, which are one nucelotide long and cannot code for
              amino acid; numeric code is 0 <br />
              <b>D:</b> Exons, which can be part of both CDS and UTR varying
              from transcript to transcript; numeric code is -2 when it is a
              part of UTR and -1 to n when it is part of CDS (same as above U
              and T cases) <br />
              <b>R:</b> Special cases of intron retention (described in the end)
            </li>
            <br />
            <li>
              <h5 style={{ color: "MediumVioletRed" }}>
                Level 2 -> To describe the inclusion frequency of exon in the
                transcripts
              </h5>
              <br /> <b>G</b>: Constitutive (or Global) exons are present in all
              the transcripts <br />
              <b>A</b>: Alternate exons are not present in all transcripts
              <br /> <b>F</b>: Exons, which when considered with their alternate
              splicing sites, can be found in every transcript <br />Numeric code here
              gives a serial identifier (from 1 to n) to the exon within each of
              the above-described categories
            </li>
            <br />
            <li>
              <h5 style={{ color: "Teal" }}>
                Level 3 -> To describe the shifting in splicing sites of exons
              </h5>
              <br />
              All these modifications are assigned names by comparing them to
              the exons of principal isoform (PI), which is chosen as the
              isoform having maximum number of coding exons (preferably reviewed) in a particular gene.{" "}
              <br />
              <b>n: </b> 3’ splice site of upstream intron is changed, leading
              to extension or constriction of exon length from its 5’ end <br />
              <b>c:</b> 5’ splice site of downstream intron is changed, leading
              to extension or constriction of exon length from its 3’ end <br />
              <b>b</b>: 5’ splice site of downstream intron and 3’ splice site
              of upstream intron both are changed simultaneously<br />Numeric code
              here also gives a serial identifier (from 0 to n) to the exon
              within each of the above described categories. <br />
              <b>0</b>: denotes the original exon as in the PI of that gene.
              Numeric code will be 0 in this case
            </li>
            <br />
            <li>
              <h5 style={{ color: "Crimson" }}>Intron Retention Cases</h5> here
              the exon identifier contains 5 parts separated by 4 colons(:)
              between them. First we have R to denote the retention event
              followed by a colon(:) and the numeric code (-2 to n) which is
              similar to that of the numeric codes in level 1, again followed by
              a colon(:) after which it contains the identifier of the exon from
              which the intron is retained followed by another colon(:) and a
              numeric code specifying the serial number(0 to n) of the retention
              event from that exon, followed by a colon(:) and lastly the
              identifier of exon upto which the intron has been retained.
            </li>
          </ul>
          <p>
            <b>Glossary:</b> Untranslated Region(UTR), Coding
            Sequence(CDS), Coding Genomic Coordinates(CGC), Frame Shift
            Events(FSE),  Alternative Transcription
            Initiation(ATI), Alternative Transcription Termination(ATT),
            Principal Isoform(PI)
          </p>
        </div>
      </div>
    );
  }
}

export default Nomenclature;
