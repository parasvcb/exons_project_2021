//export const MDgeneCard = (gene) =>
export const showalert = (state) => {
    //console.log('**',state)
    if (state.exoninf === false) {
      return (
        <div
          className="alert-light align-middle bg-white px-3 py-3"
          style={{ height: "80px" }}
          role="alert"
        >
          <strong>Click on "Show exons nomenclature" button</strong>
        </div>
      );
    }
    if (state.alert) {
      return (
        <div
          className="alert alert-light"
          role="alert"
          style={{ maxHeight: "330px", minHeight: "330px", overflow: "auto"}}
        >
          {/* {this.decodeExon(this.state.alert)} */}
          {state.alert}
        </div>
      );
    } else {
      if (state.scrollable === "fasta") {
        return (
          <div
            className="align-middle p-5"
            role="alert"
            style={{ height: "180px", alignItems: "middle", background:"WHITE" }}
          >
            <span className="align-middle">
              Change view for exon properties
            </span>
          </div>
        );
      }
      return (
        <div
          className="alert alert-light"
          role="alert"
          style={{ maxHeight: "330px", minHeight: "330px", overflow: "auto" }}
        >
          <span className="align-middle">Hover over exon</span>
        </div>
      );
    }
    //alert({ id });
  }
  export const showInfFunc=(val,props) =>{
    let pfamurl = "https://pfam.xfam.org/family/";
    if (val === "ss") {
      return (
        <div className="container">
          <h5 style={{ color: "LightSlateGray" }}>
            Predicted secondary structure of exons
          </h5>
          <div className="card-body p-1 m-0">
            <span>
              <span
                className="badge badge-pill mr-3"
                style={{ background: "purple", color: "purple", width: "50px" }}
              >
                H
              </span>
            </span>

            <span>HELIX</span>
            <br />
            <span>
              <span
                className="badge badge-pill mr-3"
                style={{ background: "yellow", color: "yellow", width: "50px" }}
              >
                E
              </span>
            </span>
            <span>SHEET</span>
            <br />
            <span>
              <span
                className="badge badge-pill mr-3"
                style={{ background: "linen", color: "black", width: "50px" }}
              >
                ------
              </span>
            </span>
            <span>COIL</span>
          </div>
        </div>
      );
    } else {
      if (val === "aaseq") {
        return (
          <div className="container">
            {" "}
            <h5 style={{ color: "LightSlateGray" }}>
              Amino acid context of exons
            </h5>
          </div>
        );
      } else {
        if (val === "dom") {
          var dom = props.domset;
          return (
            <>
              <div className="container" style={{ padding: "0px" }}>
                <h3>Domains</h3>
                <table className="table table-sm">
                  <thead>
                    <tr>
                      <th>Color</th>
                      <th>Domain Name</th>
                      <th>Pfam ID</th>
                    </tr>
                  </thead>
                  <tbody>
                    {dom.map(ds => (
                      <tr key={ds.code}>
                        <td>
                          <span
                            className="badge badge-light"
                            style={{ background: ds.color, color: ds.color }}
                          >
                            ds.code
                          </span>
                        </td>
                        <td>{ds.dId}</td>
                        <td>
                          <a target="_blank" href={pfamurl + ds.name}>
                            {ds.name}
                          </a>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </>
          );
        } else {
          if (val === "dis") {
            return (
              <div className="container">
                <h5 style={{ color: "LightSlateGray" }}>
                  Predicted Disordered Regions of exons
                </h5>
                <span
                  className="badge badge-pill mr-3"
                  style={{ background: "olive", color: "olive", width: "50px" }}
                >
                  S
                </span>&nbsp;&nbsp;
                <span>DISORDERED REGION</span>
                <br />
                <span
                  className="badge badge-pill mr-3"
                  style={{ background: "white", color: "black", width: "50px" }}
                >
                  -------
                </span>&nbsp;&nbsp;
                <span>STRUCTURED REGION</span>
              </div>
            );
          } else {
            return (
              <div className="container">
                <h5 style={{ color: "LightSlateGray" }}>Fasta view</h5>
              </div>
            );
          }
        }
      }
    }
  }