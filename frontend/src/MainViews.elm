module MainViews exposing (..)

import Array
import Dict
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events as Events exposing (onClick)
import MainHelpers exposing (queryString)
import MainTypes exposing (..)
import SearchBarTypes exposing (SearchResult, SearchResults)
import Table
import Url.Builder as UBuilder


selectClickedResult : SearchResult -> List (Html.Attribute Msg)
selectClickedResult result =
    [ class
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
            , Table.stringColumn "Species" (get "species")
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


tableCustomizations : Table.Customizations SearchResult Msg
tableCustomizations =
    let
        default =
            Table.defaultCustomizations
    in
    { default
        | tableAttrs = [ class "table table-hover table-sm table-bordered table-responsive-md" ]
        , rowAttrs = selectClickedResult
    }


viewSearchResults : Model -> List (Html Msg)
viewSearchResults { searchResultRows, resultsTableState, resultsTableQuery, searchHits } =
    [ div [ class "d-flex bg-light text-primary" ]
        [ text "Hits: ", text (Maybe.withDefault "" (Maybe.map String.fromInt searchHits)) ]
    , Table.view tableConfig resultsTableState (Maybe.withDefault [] (Maybe.map Array.toList searchResultRows))
    , a
        [ class "btn btn-outline-primary"
        , href <| "/download/" ++ (queryString searchResultRows)
        , attribute "download" "data.zip"
        ]
        [ text "download" ]
    ]
