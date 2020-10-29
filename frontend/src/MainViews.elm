module MainViews exposing (..)

import Array
import Dict
import Html exposing (..)
import Html.Attributes as Attr
import Html.Events as Events
import MainTypes exposing (..)
import SearchBarTypes exposing (SearchResult, SearchResults)
import Table


selectClickedResult : SearchResult -> List (Html.Attribute Msg)
selectClickedResult result =
    [ Attr.class
        (if result.selected then
            "table-primary"

         else
            ""
        )
    , Events.onClick (ResultClicked result)
    ]


get : String -> SearchResult -> String
get key searchResult =
    Maybe.withDefault "" (Dict.get key searchResult.data)


tableConfig : Table.Config SearchResult Msg
tableConfig =
    let
        getId =
            .id >> String.fromInt
    in
    Table.customConfig
        { toId = getId
        , toMsg = SetResultsTableState
        , columns =
            [ Table.stringColumn "Row" getId
            , Table.stringColumn "SRA Run" (get "SRR_accession")
            , Table.stringColumn "QC summary" (get "QC_summary")
            , Table.stringColumn "SRA experiment" (get "SRX_accession")
            , Table.stringColumn "SRA sample" (get "SRS_accession")
            , Table.stringColumn "SRA project" (get "SRP_accession")
            , Table.stringColumn "Sample" (get "Sample_name")
            , Table.stringColumn "Experiment" (get "GEO_series")
            ]
         , customizations = tableCustomizations
        }

tableCustomizations: Table.Customizations SearchResult Msg
tableCustomizations =
    let
        default = Table.defaultCustomizations
    in
       {default | tableAttrs = [ Attr.class "table table-hover table-sm table-bordered table-responsive" ]
       , rowAttrs = selectClickedResult
       }


--# This is the header of the current search result page
--# ['SRA run accession', 'QC summary alttext ', 'SRA experiment accession', 'SRA sample accession',
--# 'SRA project accession', 'Sample Name / GEO sample accession', 'GEO series accession', 'Experiment name']

viewSearchResults : Model -> Html Msg
viewSearchResults { searchResultRows, resultsTableState, resultsTableQuery } =
    Table.view tableConfig resultsTableState (Maybe.withDefault [] (Maybe.map Array.toList searchResultRows))
