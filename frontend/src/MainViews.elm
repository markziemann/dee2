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


selectClickedResult : SearchResult -> List (Html.Attribute Msg)
selectClickedResult result =
    [ class
        (if result.selected then
            "table-primary"

         else
            ""
        )
    , style "cursor" "pointer" -- First place not using bootstrap for style?
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
            --GEO_series
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


hideWhenTrue : String -> Bool -> String
hideWhenTrue classString true =
    if true then
        String.join " " [ classString, "invisible" ]

    else
        classString


noResultsSelected maybeRows =
    List.all identity <| Array.toList <| Array.map (.selected >> not) <| Maybe.withDefault Array.empty maybeRows


buttonOrSpinner : Bool -> Maybe (Array.Array SearchResult) -> Html Msg
buttonOrSpinner downloading rows =
    if not downloading then
        a
            [ hideWhenTrue "btn btn-outline-primary btn-block" (noResultsSelected rows) |> class
            , href <| "/download/" ++ queryString rows
            , attribute "download" "data.zip"
            , onClick DownloadRequested
            ]
            [ text "download" ]

    else
        button
            [ class "btn btn-primary btn-block"
            , attribute "type" "button"
            , attribute "disabled" "disabled"
            ]
            [ span
                [ class "spinner-border spinner-border-sm mr-2"
                , attribute "role" "status"
                , attribute "aria-hidden" "true"
                ]
                []
            , text "Loading..."
            ]


viewSearchResults : Model -> List (Html Msg)
viewSearchResults ({ searchResultRows, resultsTableState, resultsTableQuery, searchHits } as model) =
    [ div [ class "d-flex bg-light text-primary" ]
        [ text "Hits: ", text (Maybe.withDefault "" (Maybe.map String.fromInt searchHits)) ]
    , Table.view tableConfig resultsTableState (Maybe.withDefault [] (Maybe.map Array.toList searchResultRows))
    , div [class "btn-group"] [buttonOrSpinner model.downloading searchResultRows]
    ]
