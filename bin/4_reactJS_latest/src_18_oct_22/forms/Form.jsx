import React from "react";

export default ({ handleChange, handleSubmit, post }) => {
  return (
    
      // <div className="p-5" style={{background: 'linear-gradient(to right, #E9E9E9, #E8E8E8'}}>
      // <div className="p-5" style={{ borderRadius: "20% 60% 50% 30% / 90% 70% 50% 30%", background: 'linear-gradient(to right, #E9E9E9, #E8E8E8'}}>
      <div className="p-5" style={{ 
        borderRadius: "20% 60% 50% 30% / 90% 70% 50% 30%", 
        background: 'linear-gradient(to right, #f8eee7, #f4decb',
        background: '#f8eee7',
        }}>
        
        <div>
        {/* will write conetnt here */}
        <br></br>
          <h4 style={{fontFamily:"Baloo2", fontWeight: "thin", color:"#52658f"}}>
            To view the exon nomenclature and their visualization in protein coding transcripts, you need to provide gene name, or NCBI gene id (Entrez)</h4>
        </div>
        <br></br>
        <br></br>
        <div
          className="form-container position-relative" align="left"
          style = {{
            display: "flex",
          }}
        >

          <form>
            <div className="mb-3">
            <label for="exampleInputEmail1" className="form-label" style={{color:"#333a56"}}>Please Enter here</label>
              <input
                className="col-sm-12"
                name="name"
                onChange={handleChange}
                type="text"
                value={post.name}
                // placeholder="Gene Name or NCBI gene Id"
                // width="100px"
                // height="30px"
              />
              <div id="passwordHelpBlock" class="form-text">
                  Your input should be either numeric (Entrez identifier) or name of gene.
              </div>
            </div>
            <br></br>
        
            <button className="btn" onClick={handleSubmit} style={{background:"#DCAE1d",color:"White"}}>
              {/* style={{background:"#88D317",color:"White"}} */}
              Submit
            </button>
          </form>
        </div>
      </div>   
  );
};
