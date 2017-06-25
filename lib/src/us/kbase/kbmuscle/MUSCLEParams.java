
package us.kbase.kbmuscle;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: MUSCLE_Params</p>
 * <pre>
 * MUSCLE Input Params
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace_name",
    "desc",
    "input_ref",
    "output_name",
    "maxiters",
    "maxhours"
})
public class MUSCLEParams {

    @JsonProperty("workspace_name")
    private String workspaceName;
    @JsonProperty("desc")
    private String desc;
    @JsonProperty("input_ref")
    private String inputRef;
    @JsonProperty("output_name")
    private String outputName;
    @JsonProperty("maxiters")
    private Long maxiters;
    @JsonProperty("maxhours")
    private Double maxhours;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace_name")
    public String getWorkspaceName() {
        return workspaceName;
    }

    @JsonProperty("workspace_name")
    public void setWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
    }

    public MUSCLEParams withWorkspaceName(String workspaceName) {
        this.workspaceName = workspaceName;
        return this;
    }

    @JsonProperty("desc")
    public String getDesc() {
        return desc;
    }

    @JsonProperty("desc")
    public void setDesc(String desc) {
        this.desc = desc;
    }

    public MUSCLEParams withDesc(String desc) {
        this.desc = desc;
        return this;
    }

    @JsonProperty("input_ref")
    public String getInputRef() {
        return inputRef;
    }

    @JsonProperty("input_ref")
    public void setInputRef(String inputRef) {
        this.inputRef = inputRef;
    }

    public MUSCLEParams withInputRef(String inputRef) {
        this.inputRef = inputRef;
        return this;
    }

    @JsonProperty("output_name")
    public String getOutputName() {
        return outputName;
    }

    @JsonProperty("output_name")
    public void setOutputName(String outputName) {
        this.outputName = outputName;
    }

    public MUSCLEParams withOutputName(String outputName) {
        this.outputName = outputName;
        return this;
    }

    @JsonProperty("maxiters")
    public Long getMaxiters() {
        return maxiters;
    }

    @JsonProperty("maxiters")
    public void setMaxiters(Long maxiters) {
        this.maxiters = maxiters;
    }

    public MUSCLEParams withMaxiters(Long maxiters) {
        this.maxiters = maxiters;
        return this;
    }

    @JsonProperty("maxhours")
    public Double getMaxhours() {
        return maxhours;
    }

    @JsonProperty("maxhours")
    public void setMaxhours(Double maxhours) {
        this.maxhours = maxhours;
    }

    public MUSCLEParams withMaxhours(Double maxhours) {
        this.maxhours = maxhours;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((("MUSCLEParams"+" [workspaceName=")+ workspaceName)+", desc=")+ desc)+", inputRef=")+ inputRef)+", outputName=")+ outputName)+", maxiters=")+ maxiters)+", maxhours=")+ maxhours)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
